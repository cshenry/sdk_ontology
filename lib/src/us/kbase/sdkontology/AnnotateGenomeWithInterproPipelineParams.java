
package us.kbase.sdkontology;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: AnnotateGenomeWithInterproPipelineParams</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace",
    "genome_workspace",
    "genome_id",
    "genome_output_id"
})
public class AnnotateGenomeWithInterproPipelineParams {

    @JsonProperty("workspace")
    private String workspace;
    @JsonProperty("genome_workspace")
    private String genomeWorkspace;
    @JsonProperty("genome_id")
    private String genomeId;
    @JsonProperty("genome_output_id")
    private String genomeOutputId;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace")
    public String getWorkspace() {
        return workspace;
    }

    @JsonProperty("workspace")
    public void setWorkspace(String workspace) {
        this.workspace = workspace;
    }

    public AnnotateGenomeWithInterproPipelineParams withWorkspace(String workspace) {
        this.workspace = workspace;
        return this;
    }

    @JsonProperty("genome_workspace")
    public String getGenomeWorkspace() {
        return genomeWorkspace;
    }

    @JsonProperty("genome_workspace")
    public void setGenomeWorkspace(String genomeWorkspace) {
        this.genomeWorkspace = genomeWorkspace;
    }

    public AnnotateGenomeWithInterproPipelineParams withGenomeWorkspace(String genomeWorkspace) {
        this.genomeWorkspace = genomeWorkspace;
        return this;
    }

    @JsonProperty("genome_id")
    public String getGenomeId() {
        return genomeId;
    }

    @JsonProperty("genome_id")
    public void setGenomeId(String genomeId) {
        this.genomeId = genomeId;
    }

    public AnnotateGenomeWithInterproPipelineParams withGenomeId(String genomeId) {
        this.genomeId = genomeId;
        return this;
    }

    @JsonProperty("genome_output_id")
    public String getGenomeOutputId() {
        return genomeOutputId;
    }

    @JsonProperty("genome_output_id")
    public void setGenomeOutputId(String genomeOutputId) {
        this.genomeOutputId = genomeOutputId;
    }

    public AnnotateGenomeWithInterproPipelineParams withGenomeOutputId(String genomeOutputId) {
        this.genomeOutputId = genomeOutputId;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((("AnnotateGenomeWithInterproPipelineParams"+" [workspace=")+ workspace)+", genomeWorkspace=")+ genomeWorkspace)+", genomeId=")+ genomeId)+", genomeOutputId=")+ genomeOutputId)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
