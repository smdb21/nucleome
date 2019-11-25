package edu.scripps.yates.nucleome.turboID;

import gnu.trove.map.hash.TObjectDoubleHashMap;

public class ProteinFromIP {
	private final TObjectDoubleHashMap<IPExperimentType> normalizedIntensities = new TObjectDoubleHashMap<IPExperimentType>();
	private final TObjectDoubleHashMap<IPExperimentType> originalIntensities = new TObjectDoubleHashMap<IPExperimentType>();
	private final String acc;
	private final boolean isTransmembrane;
	private int spc;
	private final String gene;

	public ProteinFromIP(String acc, String gene, boolean isTransmembrane) {
		this.acc = acc;
		this.isTransmembrane = isTransmembrane;
		this.gene = gene;

	}

	public boolean addNormalizedIntensity(double intensity, IP_Channel_Norm channel) {
		return addNormalizedIntensity(intensity, channel.getExpType());
	}

	public boolean addNormalizedIntensity(double intensity, IPExperimentType expType) {
		if (Double.isNaN(intensity)) {
			return false;
		}
		normalizedIntensities.put(expType, intensity);

		return true;
	}

	public String getAcc() {
		return acc;
	}

	public boolean isTransmembrane() {
		return isTransmembrane;
	}

	public Double getNormalizedIntensity(IPExperimentType experiment) {
		return normalizedIntensities.get(experiment);
	}

	public Double getOriginalIntensity(IPExperimentType experiment) {
		return originalIntensities.get(experiment);
	}

	public Double getNormalizedIntensity(IP_Channel_Norm channel) {
		return normalizedIntensities.get(channel.getExpType());
	}

	public Double getOriginalIntensity(IP_Channel_Ori channel) {
		return originalIntensities.get(channel.getExpType());
	}

	public String getGene() {
		return this.gene;
	}

	public boolean addOriginalIntensity(Double intensity, IP_Channel_Ori channel) {
		return addOriginalIntensity(intensity, channel.getExpType());

	}

	public boolean addOriginalIntensity(Double intensity, IPExperimentType expType) {
		if (Double.isNaN(intensity)) {
			return false;
		}

		originalIntensities.put(expType, intensity);

		return true;

	}

	public int getSpc() {
		return spc;
	}

	public void setSpc(int spc) {
		this.spc = spc;
	}
}
