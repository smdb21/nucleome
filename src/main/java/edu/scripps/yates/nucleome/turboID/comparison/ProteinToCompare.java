package edu.scripps.yates.nucleome.turboID.comparison;

import edu.scripps.yates.nucleome.turboID.TurboIDExperimentType;
import gnu.trove.map.hash.TObjectDoubleHashMap;

public class ProteinToCompare {
	private final String acc;
	private final String gene;
	private final String description;
	private final Integer clusterNum;
	private final double pvalue;
	private final TObjectDoubleHashMap<TurboIDExperimentType> quantMap = new TObjectDoubleHashMap<TurboIDExperimentType>();

	public ProteinToCompare(String acc, String gene, String description, Integer clusterNum, double pvalue) {
		this.acc = acc;
		this.gene = gene;
		this.description = description;
		this.clusterNum = clusterNum;
		this.pvalue = pvalue;
	}

	public void addQuant(TurboIDExperimentType bait, double quant) {
		quantMap.put(bait, quant);
	}

	public double getQuant(TurboIDExperimentType bait) {
		return quantMap.get(bait);
	}

	public String getAcc() {
		return acc;
	}

	public String getGene() {
		return gene;
	}

	public String getDescription() {
		return description;
	}

	public Integer getClusterNum() {
		return clusterNum;
	}

	public double getPvalue() {
		return pvalue;
	}
}
