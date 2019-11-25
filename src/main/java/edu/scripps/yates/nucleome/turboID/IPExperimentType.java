package edu.scripps.yates.nucleome.turboID;

public enum IPExperimentType {
	NU, BKGD300, EMD150a(true), EMD300a(true), EMD150b(true), EMD300b(true), Bt300;
	private final boolean isBait;

	private IPExperimentType() {
		this(false);
	}

	private IPExperimentType(boolean isBait) {
		this.isBait = isBait;
	}

	public boolean isBait() {
		return this.isBait;
	}
}
