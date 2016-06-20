package sdk_ontology::sdk_ontologyImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = "0.1.0";

=head1 NAME

sdk_ontology

=head1 DESCRIPTION

A KBase module: sdk_ontology
This module convert given KBase annotations of a genome to GO terms.

=cut

#BEGIN_HEADER
use Bio::KBase::AuthToken;
use Bio::KBase::workspace::Client;
use Config::IniFiles;
use Data::Dumper;
use File::Path;
use DateTime;
use Cwd;
our $EC_PATTERN = qr/\(\s*E\.?C\.?(?:\s+|:)(\d\.(?:\d+|-)\.(?:\d+|-)\.(?:n?\d+|-)\s*)\)/;
our $currentcontext;

#Initialization function for call
sub util_initialize_call {
	my ($self,$params,$ctx) = @_;
	print("Starting ".$ctx->method()." method.\n");
	$currentcontext = $ctx;
	return $params;
}

sub util_currentuser {
	return $currentcontext->user_id();
}

sub util_token {
	return $currentcontext->token();
}

sub util_provenance {
	return $currentcontext->provenance();
}

sub util_error {
	my ($self,$message) = @_;
	Carp::confess($message);
}

sub util_to_fasta {
	my ($self,$seqName, $seq, $len) = @_;
	# default to 80 characters of sequence per line
	$len = 80 unless $len;
	my $formatted_seq = ">$seqName\n";
	while (my $chunk = substr($seq, 0, $len, "")) {
		$formatted_seq .= "$chunk\n";
	}
	return $formatted_seq;
}

sub util_scratchdir {
	my ($self) = @_;
	return $self->{_scratchdir};
}

sub util_validate_args {
	my ($self,$args,$mandatoryArguments,$optionalArguments,$substitutions) = @_;
	print "Retrieving input parameters.\n";
	if (!defined($args)) {
	    $args = {};
	}
	if (ref($args) ne "HASH") {
		$self->util_error("Arguments not hash");	
	}
	if (defined($substitutions) && ref($substitutions) eq "HASH") {
		foreach my $original (keys(%{$substitutions})) {
			$args->{$original} = $args->{$substitutions->{$original}};
		}
	}
	if (defined($mandatoryArguments)) {
		for (my $i=0; $i < @{$mandatoryArguments}; $i++) {
			if (!defined($args->{$mandatoryArguments->[$i]})) {
				push(@{$args->{_error}},$mandatoryArguments->[$i]);
			}
		}
	}
	$self->util_error("Mandatory arguments ".join("; ",@{$args->{_error}})." missing.") if (defined($args->{_error}));
	if (defined($optionalArguments)) {
		foreach my $argument (keys(%{$optionalArguments})) {
			if (!defined($args->{$argument})) {
				$args->{$argument} = $optionalArguments->{$argument};
			}
		}	
	}
	return $args;
}

sub searchname
{
    my $sn = $_[0];
    $sn =~ s/^\s+//;
    $sn =~ s/_//g;
    $sn =~ s/-//g;
    $sn =~ s/,//g;
    $sn =~ tr/A-Z/a-z/;
    $sn =~ s/[\s]//g;
    $sn =~ s/(?<=\(ec)[^)]+[^(]+(?=\))//g;
    $sn =~ s/(?<=\(tc)[^)]+[^(]+(?=\))//g;
    return $sn;
}

sub splitFunc
{
    my $fn = $_[0];
    my @splitFunc;
    if ($fn =~ /\// || $fn =~ /;/ || $fn =~ /@/){
        @splitFunc = split /[;@\/]+/, $fn;
        return \@splitFunc;
    }
    else{
        push (@splitFunc, $fn);
        return \@splitFunc;
    }
}

sub searchec
{
    my $sn = $_[0];
    $sn =~ s/://g;
    $sn =~ tr/A-Z/a-z/;
    $sn =~ s/[\s]//g;
    $sn =~ s/^\s+//;
    return $sn;
}

sub featureTranslate{
    my ($genome, $ontTr, $ontRef, $ont_tr, $clear) = @_;
    my $func_list = $genome->{features};
    my %selectedRoles;
    my %termName;
    my %termId;

    foreach my $k (keys $ontTr){
        my $r = $ontTr->{$k}->{name};
        my $eq = $ontTr->{$k}->{equiv_terms};
        my $mRole = searchname ($r);
        my @ecArr;
        if ($ont_tr eq "ec2go"){
            my $ecsn = searchec($k);
            $mRole=$ecsn;
        }
            my @tempMR;
            for (my $i=0; $i<@$eq; $i++){
                my $e_name = $eq->[$i]->{equiv_name};
                my $e_term = $eq->[$i]->{equiv_term};
                $termName{$e_name} = $e_term;
                $termId{$mRole} = [$r,$k];
                push (@tempMR, $e_name);
            }
            $selectedRoles{$mRole} = \@tempMR;
    }
    print "Following feature annotations were translated in the genome\n";
    my $local_time = localtime ();
    my $vs = version ();
    my $changeRoles =0;
    for (my $j =0; $j< @$func_list; $j++){

        if ($clear == 1){
            $func_list->[$j]->{ontology_terms} = {};

        }
        my $func = $func_list->[$j]->{function};
        my $funcId = $func_list->[$j]->{id};
        my $splitArr = splitFunc($func);
        my $count_flag =0;
        foreach my $fr (@$splitArr){
            $count_flag++;
            my $sn;
            if ($ont_tr eq "ec2go"){
                my $ecNum;
                if ($fr =~ /(.+?)\s*$EC_PATTERN\s*(.*)/) {
                    $ecNum = $2;
                    $sn = searchec ("EC:$ecNum");
                }
            }
            else{
                 $sn = searchname($fr);
            }
            ########## for unitprot and ec partial mappings being considered######
            if ($ont_tr eq "uniprotkb_kw2go" || $ont_tr eq "ec2go" ){
                while (my($key, $value) = each %selectedRoles) {
                    if (-1 != index($sn, $key)) {
                        $sn =$key;

                        if ( exists $selectedRoles{$sn}  && !defined ($func_list->[$j]->{ontology_terms}->{GO}) ){
                            my $nrL = $selectedRoles{$sn};
                            my @tempA;
                            for (my $i=0; $i< @$nrL; $i++){
                                my $ontEvi = {
                                    method => "Remap annotations based on Ontology translation table",
                                    method_version => $vs,
                                    timestamp => $local_time,
                                    translation_provenance => [],
                                    alignment_evidence => []
                                };
                                my $ontData ={
                                    id => "",
                                    ontology_ref =>"",
                                    term_lineage => [],
                                    term_name => "",
                                    evidence => []
                                };
                                push(@tempA, $nrL->[$i]);
                                $ontData->{id} = $termName{$nrL->[$i]}; #$termName{$nrL->[$i]};
                                $ontData->{ontology_ref} = "dictionary ref";
                                $ontData->{term_name} = $nrL->[$i];#$nrL->[$i];
                                $ontEvi->{translation_provenance} = [$ontRef, $termId{$sn}->[1], $termId{$sn}->[0]];
                                push (@{$ontData->{evidence}},$ontEvi);
                                $func_list->[$j]->{ontology_terms}->{GO}->{$termName{$nrL->[$i]}} = $ontData;
                            }
                            my $joinStr = join (" | ", @tempA);
                            print "$funcId\t$fr\t-\t$joinStr\n";
                            #print &Dumper ($func_list->[$j]->{ontology_terms});
                            $changeRoles++;
                        }
                        elsif ( exists $selectedRoles{$sn} && defined ($func_list->[$j]->{ontology_terms}->{GO}) && $count_flag <= 1 ) {
                                my $new_term = $func_list->[$j]->{ontology_terms}->{GO};

                                my $nrL = $selectedRoles{$sn};
                                my @tempA;
                                for (my $i=0; $i< @$nrL; $i++){
                                    my $ontEvi = {
                                        method => "Remap annotations based on Ontology translation table",
                                        method_version => $vs,
                                        timestamp => $local_time,
                                        translation_provenance => [],
                                        alignment_evidence => []
                                    };
                                    my $ontData ={
                                        id => "",
                                        ontology_ref =>"",
                                        term_lineage => [],
                                        term_name => "",
                                        evidence => []
                                    };
                                    push(@tempA, $nrL->[$i]);

                                    if (exists $new_term->{$termName{$nrL->[$i]}}){
                                        $ontEvi->{translation_provenance} = [$ontRef, $termId{$sn}->[1], $termId{$sn}->[0]];
                                        push (@{$new_term->{$termName{$nrL->[$i]}}->{evidence}}, $ontEvi);
                                    }
                                    else{
                                        $ontData->{id} = $termName{$nrL->[$i]}; #$termName{$nrL->[$i]};
                                        $ontData->{ontology_ref} = "dictionary ref";
                                        $ontData->{term_name} = $nrL->[$i];#$nrL->[$i];
                                        $ontEvi->{translation_provenance} = [$ontRef, $termId{$sn}->[1], $termId{$sn}->[0]];
                                        push (@{$ontData->{evidence}},$ontEvi);
                                        $func_list->[$j]->{ontology_terms}->{GO}->{$termName{$nrL->[$i]}} = $ontData;
                                    }
                                }
                                my $joinStr = join (" | ", @tempA);
                                print "$funcId\t$fr\t-\t$joinStr\n";
                                $changeRoles++;
                                #print &Dumper ( $func_list->[$j]->{ontology_terms});
                                #die;
                        }
                        else{
                            next;
                        }
#################################end
                    }
                }
            }
            else{

                if ( exists $selectedRoles{$sn}  && !defined ($func_list->[$j]->{ontology_terms}->{GO}) ){
                    my $nrL = $selectedRoles{$sn};
                        my @tempA;
                        for (my $i=0; $i< @$nrL; $i++){
                            my $ontEvi = {
                                method => "Remap annotations based on Ontology translation table",
                                method_version => $vs,
                                timestamp => $local_time,
                                translation_provenance => [],
                                alignment_evidence => []
                            };
                            my $ontData ={
                                id => "",
                                ontology_ref =>"",
                                term_lineage => [],
                                term_name => "",
                                evidence => []
                            };
                            push(@tempA, $nrL->[$i]);
                            $ontData->{id} = $termName{$nrL->[$i]}; #$termName{$nrL->[$i]};
                            $ontData->{ontology_ref} = "dictionary ref";
                            $ontData->{term_name} = $nrL->[$i];#$nrL->[$i];
                            $ontEvi->{translation_provenance} = [$ontRef, $termId{$sn}->[1], $termId{$sn}->[0]];
                            push (@{$ontData->{evidence}},$ontEvi);
                            $func_list->[$j]->{ontology_terms}->{GO}->{$termName{$nrL->[$i]}} = $ontData;
                        }
                        my $joinStr = join (" | ", @tempA);
                        print "$funcId\t$fr\t-\t$joinStr\n";
                        #print &Dumper ($func_list->[$j]->{ontology_terms});
                        $changeRoles++;
                }
                elsif ( exists $selectedRoles{$sn} && defined ($func_list->[$j]->{ontology_terms}->{GO}) && $count_flag <= 1 ) {
                        my $new_term = $func_list->[$j]->{ontology_terms};
                        my $nrL = $selectedRoles{$sn};
                        my @tempA;
                        for (my $i=0; $i< @$nrL; $i++){
                            my $ontEvi = {
                                method => "Remap annotations based on Ontology translation table",
                                method_version => $vs,
                                timestamp => $local_time,
                                translation_provenance => [],
                                alignment_evidence => []
                            };
                            my $ontData ={
                                id => "",
                                ontology_ref =>"",
                                term_lineage => [],
                                term_name => "",
                                evidence => []
                            };
                            push(@tempA, $nrL->[$i]);

                            if (exists $new_term->{$termName{$nrL->[$i]}}){
                                $ontEvi->{translation_provenance} = [$ontRef, $termId{$sn}->[1], $termId{$sn}->[0]];
                                push (@{$new_term->{$termName{$nrL->[$i]}}->{evidence}}, $ontEvi);
                            }
                            else{
                                $ontData->{id} = $termName{$nrL->[$i]}; #$termName{$nrL->[$i]};
                                $ontData->{ontology_ref} = "dictionary ref";
                                $ontData->{term_name} = $nrL->[$i];#$nrL->[$i];
                                $ontEvi->{translation_provenance} = [$ontRef, $termId{$sn}->[1], $termId{$sn}->[0]];
                                push (@{$ontData->{evidence}},$ontEvi);
                                $func_list->[$j]->{ontology_terms}->{GO}->{$termName{$nrL->[$i]}} = $ontData;
                            }
                        }
                        my $joinStr = join (" | ", @tempA);
                        print "$funcId\t$fr\t-\t$joinStr\n";
                        $changeRoles++;
                        #print &Dumper ( $func_list->[$j]->{ontology_terms});
                        #die;
                }
                else{
                    next;
                }
            }
        } #foreach
    }#for
    print "\nTotal of $changeRoles feature annotations were translated \n";
}

sub func_annotate_genome_with_interpro_pipeline {
	my ($self,$params) = @_;
    $params = $self->util_validate_args($params,["workspace","genome_id","genome_output_id"],{
    	genome_workspace => $params->{workspace},
    });
    my $annofunc = "Annotate Genome with InterPro Pipeline";
  	my $timestamp = DateTime->now()->datetime();
    #Step 1: Get genome from workspace
    my $wsClient = Bio::KBase::workspace::Client->new($self->{'workspace-url'},token=>$self->util_token());
    my $genome = $wsClient->get_objects([{workspace=>$params->{genome_workspace},name=>$params->{genome_id}}])->[0]{data};
    #Step 2: Print protein FASTA file
    File::Path::mkpath $self->util_scratchdir()
    my $filename = $self->util_scratchdir()."/protein.fa";
    open ( my $fh, ">", $filename) || $self->util_error("Failure to open file: $filename, $!");
    my $genehash = {};
    foreach my $gene (@{$genome->{features}}) {
    	if (defined($gene->{protein_translation}) && length($gene->{protein_translation}) > 0) {
    		$genehash->{$gene->{id}} = $gene;
    		print $fh $self->util_to_fasta($gene->{id}, $gene->{protein_translation});
    	}
    }
    close($fh);
    #Step 3: Run interpro
    my $orig_cwd = cwd;
    chdir $self->util_scratchdir();
    system("interproscan.sh -i protein.fa -f tsv -o protein.tsv --disable-precalc -iprscan -iprlookup -hm");
    chdir $orig_cwd;
    #Step 4: Parsing interpro results
    $filename = $self->util_scratchdir()."/protein.tsv";
    my $numftr = 0;
    my $numdomains = 0;
    my $domainhash = {};
    my $ftrhash = {};
    open ( my $fh, "<", $filename) || $self->util_error("Failure to open file: $filename, $!");
    while (my $line = <$fh>) {
    	chomp($line);
    	my $array = [split(/\t/,$line)];
    	if (@{$array} < 13 || !defined($genehash->{$array->[0]})) {
    		next;
    	}
    	if (!defined($ftrhash->{$array->[0]})) {
    		$ftrhash->{$array->[0]} = 1;
    		$numftr++;
    	}
    	if (!defined($domainhash->{$array->[11]})) {
    		$domainhash->{$array->[11]} = 1;
    		$numdomains++;
    	}
    	my $ftr = $genehash->{$array->[0]};
    	if (!defined($ftr->{ontology_terms}->{InterPro}->{$array->[11]})) {
			$ftr->{ontology_terms}->{InterPro}->{$array->[11]} = {
				 evidence => [],
				 id => $array->[11],
				 term_name => $array->[12],
				 ontology_ref => "7537/36/2",#TODO: Need to make an interpro ontology and then set this ref to that
				 term_lineage => [],
			};
		}
		my $found = 0;
		for (my $k=0; $k < @{$ftr->{ontology_terms}->{InterPro}->{$array->[11]}->{evidence}}; $k++) {
			if ($ftr->{ontology_terms}->{InterPro}->{$array->[11]}->{evidence}->[$k]->{method} eq $annofunc) {
				$ftr->{ontology_terms}->{InterPro}->{$array->[11]}->{evidence}->[$k]->{timestamp} = $timestamp;
				$ftr->{ontology_terms}->{InterPro}->{$array->[11]}->{evidence}->[$k]->{method_version} = $version;
				$ftr->{ontology_terms}->{InterPro}->{$array->[11]}->{evidence}->[$k]->{alignment_evidence} = [[$array->[6],$array->[7],abs($array->[7]-$array->[6]),$array->[8]]];
				$found = 1;
				last;
			}
		}
		if ($found == 0) {
			push(@{$ftr->{ontology_terms}->{SSO}->{$funchash->{$rolename}->{id}}->{evidence}},{
				method => $annofunc,
				method_version => $self->version,
				timestamp => $timestamp,
				$ftr->{ontology_terms}->{InterPro}->{$array->[11]}->{evidence}->[$k]->{alignment_evidence} = [[$array->[6],$array->[7],abs($array->[7]-$array->[6]),$array->[8]]];
			});
		}
    }
    close($fh);
    #Step 5: Saving the genome and report
    my $info = $wsClient->save_objects({
		'workspace'=>$params->{workspace},
		'objects'=>[{
			'type'=>'KBaseGenomes.Genome',
			'data'=>$genome,
			'name'=>$params->{genome_output_id},
			'provenance'=>$self->util_provenance()
		}]
	});
    my $reportObj = {
		'objects_created' => [$info->[6]."/".$info->[0]."/".$info->[4]],
		'text_message' => $numftr." annotated with ".$numdomains." distinct interpro domains by interpro scan!"
	};
    $info = $wsClient->save_objects({
    	workspace => $params->{workspace},
    	objects => [{
    		type => "KBaseReport.Report",
    		data => $reportObj,
    		name => $params->{genome_output_id}.".annotate_genome_with_interpro_pipeline.report",
    		hidden => 1,
    		provenance => $self->util_provenance(),
    		meta => {}
    	}]
    });
   	return {
		report_name => $params->{genome_output_id}.".annotate_genome_with_interpro_pipeline.report",
		ws_report_id => $params->{workspace}.'/'.$params->{genome_output_id}.".annotate_genome_with_interpro_pipeline.report";
	};
}

##################################

#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR

    my $config_file = $ENV{ KB_DEPLOYMENT_CONFIG };
    my $cfg = Config::IniFiles->new(-file=>$config_file);
    my $wsInstance = $cfg->val('sdk_ontology','workspace-url');
    $self->{_scratchdir} = $cfg->val('sdk_ontology','scratch');
    die "no workspace-url defined" unless $wsInstance;

    $self->{'workspace-url'} = $wsInstance;

    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}

=head1 METHODS



=head2 annotationtogo

  $output = $obj->annotationtogo($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a sdk_ontology.ElectronicAnnotationParams
$output is a sdk_ontology.ElectronicAnnotationResults
ElectronicAnnotationParams is a reference to a hash where the following keys are defined:
	workspace has a value which is a string
	input_genome has a value which is a string
	ontology_translation has a value which is a string
	translation_behavior has a value which is a string
	custom_translation has a value which is a string
	clear_existing has a value which is a string
	output_genome has a value which is a string
ElectronicAnnotationResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	output_genome_ref has a value which is a string
	n_total_features has a value which is an int
	n_features_mapped has a value which is an int

</pre>

=end html

=begin text

$params is a sdk_ontology.ElectronicAnnotationParams
$output is a sdk_ontology.ElectronicAnnotationResults
ElectronicAnnotationParams is a reference to a hash where the following keys are defined:
	workspace has a value which is a string
	input_genome has a value which is a string
	ontology_translation has a value which is a string
	translation_behavior has a value which is a string
	custom_translation has a value which is a string
	clear_existing has a value which is a string
	output_genome has a value which is a string
ElectronicAnnotationResults is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	output_genome_ref has a value which is a string
	n_total_features has a value which is an int
	n_features_mapped has a value which is an int


=end text



=item Description



=back

=cut

sub annotationtogo
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotationtogo:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotationtogo');
    }

    my $ctx = $sdk_ontology::sdk_ontologyServer::CallContext;
    my($output);
    #BEGIN annotationtogo

    print("Starting sdk_ontology method...\n\n");

    if (!exists $params->{'workspace'}) {
        die "Parameter workspace is not set in input arguments";
    }
    my $workspace_name=$params->{'workspace'};


    if (!exists $params->{'input_genome'}) {
        die "Parameter input_genome is not set in input arguments";
    }
    my $input_gen=$params->{'input_genome'};


    if (!exists $params->{'ontology_translation'}) {
        die "Parameter ontology_translation is not set in input arguments";
    }
    my $ont_tr=$params->{'ontology_translation'};

    if (!exists $params->{'translation_behavior'}) {
        die "Parameter translation_behavior is not set in input arguments";
    }
    my $trns_bh=$params->{'translation_behavior'};

    if (!exists $params->{'clear_existing'}) {
        die "Parameter clear_existing is not set in input arguments";
    }
    my $cl_ex=$params->{'clear_existing'};

    my $cus_tr;
    if (!exists $params->{'custom_translation'} && $ont_tr eq "custom") {
        die "Provide the custom translation table as an input\n\n";
    }
    elsif (exists $params->{'custom_translation'} && $ont_tr ne "custom"){
        print "Using the selected translational table from the dropdown..\n\n"
    }
    else{
        $cus_tr=$params->{'custom_translation'};
        print "Using the custom translational table..\n\n"
    }

    if (!exists $params->{'output_genome'}) {
        die "Parameter output_genome is not set in input arguments";
    }
    my $outGenome=$params->{'output_genome'};

    my $token=$ctx->token;
    my $provenance=$ctx->provenance;
    my $wsClient=Bio::KBase::workspace::Client->new($self->{'workspace-url'},token=>$token);
    my $genome=undef;
    my $ontTr=undef;
    my $cusTr=undef;
    my $ontWs = "KBaseOntology";

    if (defined $cus_tr && $ont_tr eq "custom"){
        $ontWs=$workspace_name;
        $ont_tr=$cus_tr;
    }
=head
    else{

        die "Custome translationial table is not provided\n\n";
    }
=cut
    eval {
        $genome=$wsClient->get_objects([{workspace=>$workspace_name,name=>$input_gen}])->[0]{data};
        $ontTr=$wsClient->get_objects([{workspace=>$ontWs,name=>$ont_tr}])->[0];#{data}{translation};
    };
    if ($@) {
        die "Error loading ontology translation object from workspace:\n".$@;
    }
    my $ontRef = $ontTr->{info}->[6]."/".$ontTr->{info}->[0]."/".$ontTr->{info}->[4];

    if ( ($ont_tr eq "sso2go" || $ont_tr eq "interpro2go" || $ont_tr eq "custom" || $ont_tr eq "uniprotkb_kw2go" || $ont_tr eq "ec2go" )  && ($trns_bh eq "featureOnly") ){
    featureTranslate($genome, $ontTr->{data}->{translation}, $ontRef, $ont_tr, $cl_ex);
    }
    my $obj_info_list = undef;
    eval {
        $obj_info_list = $wsClient->save_objects({
            'workspace'=>$workspace_name,
            'objects'=>[{
                'type'=>'KBaseGenomes.Genome',
                'data'=>$genome,
                'name'=>$outGenome,
                'provenance'=>$provenance
            }]
        });
    };
    if ($@) {
        die "Error saving modified genome object to workspace:\n".$@;
    }
    my $info = $obj_info_list->[0];
    print "\nMethod sucuessfully completed\n";
    print("saved:".Dumper($info)."\n");
    $output = { 'Ontology Translator' => $obj_info_list};

    #END annotationtogo
    my @_bad_returns;
    (ref($output) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"output\" (value was \"$output\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotationtogo:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotationtogo');
    }
    return($output);
}




=head2 annotate_genome_with_interpro_pipeline

  $output = $obj->annotate_genome_with_interpro_pipeline($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a sdk_ontology.AnnotateGenomeWithInterproPipelineParams
$output is a sdk_ontology.StandardFunctionOutput
AnnotateGenomeWithInterproPipelineParams is a reference to a hash where the following keys are defined:
	workspace has a value which is a sdk_ontology.workspace_name
	genome_workspace has a value which is a sdk_ontology.workspace_name
	genome_id has a value which is a sdk_ontology.genome_id
	genome_output_id has a value which is a sdk_ontology.genome_id
workspace_name is a string
genome_id is a string
StandardFunctionOutput is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a sdk_ontology.ws_report_id
ws_report_id is a string

</pre>

=end html

=begin text

$params is a sdk_ontology.AnnotateGenomeWithInterproPipelineParams
$output is a sdk_ontology.StandardFunctionOutput
AnnotateGenomeWithInterproPipelineParams is a reference to a hash where the following keys are defined:
	workspace has a value which is a sdk_ontology.workspace_name
	genome_workspace has a value which is a sdk_ontology.workspace_name
	genome_id has a value which is a sdk_ontology.genome_id
	genome_output_id has a value which is a sdk_ontology.genome_id
workspace_name is a string
genome_id is a string
StandardFunctionOutput is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a sdk_ontology.ws_report_id
ws_report_id is a string


=end text



=item Description



=back

=cut

sub annotate_genome_with_interpro_pipeline
{
    my $self = shift;
    my($params) = @_;

    my @_bad_arguments;
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_genome_with_interpro_pipeline:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_genome_with_interpro_pipeline');
    }

    my $ctx = $sdk_ontology::sdk_ontologyServer::CallContext;
    my($output);
    #BEGIN annotate_genome_with_interpro_pipeline
    $self->util_initialize_call($params,$ctx);
	$output = $self->func_annotate_genome_with_interpro_pipeline($params);
    #END annotate_genome_with_interpro_pipeline
    my @_bad_returns;
    (ref($output) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"output\" (value was \"$output\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_genome_with_interpro_pipeline:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_genome_with_interpro_pipeline');
    }
    return($output);
}




=head2 version 

  $return = $obj->version()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module version. This is a Semantic Versioning number.

=back

=cut

sub version {
    return $VERSION;
}

=head1 TYPES



=head2 bool

=over 4



=item Description

A binary boolean


=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 genome_id

=over 4



=item Description

A string representing a Genome id.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 workspace_name

=over 4



=item Description

A string representing a workspace name.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 ws_report_id

=over 4



=item Description

The workspace ID for a Report object
@id ws KBaseReport.Report


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 ElectronicAnnotationParams

=over 4



=item Description

workspace - the name of the workspace for input/output
input_genome - reference to the input genome object
ontology_translation - optional reference to user specified ontology translation map
output_genome - the name of the mapped genome annotation object

@optional ontology_translation


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace has a value which is a string
input_genome has a value which is a string
ontology_translation has a value which is a string
translation_behavior has a value which is a string
custom_translation has a value which is a string
clear_existing has a value which is a string
output_genome has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace has a value which is a string
input_genome has a value which is a string
ontology_translation has a value which is a string
translation_behavior has a value which is a string
custom_translation has a value which is a string
clear_existing has a value which is a string
output_genome has a value which is a string


=end text

=back



=head2 ElectronicAnnotationResults

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
output_genome_ref has a value which is a string
n_total_features has a value which is an int
n_features_mapped has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
output_genome_ref has a value which is a string
n_total_features has a value which is an int
n_features_mapped has a value which is an int


=end text

=back



=head2 StandardFunctionOutput

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a sdk_ontology.ws_report_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a sdk_ontology.ws_report_id


=end text

=back



=head2 AnnotateGenomeWithInterproPipelineParams

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
workspace has a value which is a sdk_ontology.workspace_name
genome_workspace has a value which is a sdk_ontology.workspace_name
genome_id has a value which is a sdk_ontology.genome_id
genome_output_id has a value which is a sdk_ontology.genome_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
workspace has a value which is a sdk_ontology.workspace_name
genome_workspace has a value which is a sdk_ontology.workspace_name
genome_id has a value which is a sdk_ontology.genome_id
genome_output_id has a value which is a sdk_ontology.genome_id


=end text

=back



=cut

1;
