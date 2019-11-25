package edu.scripps.yates.nucleome.turboID;

import java.util.ArrayList;
import java.util.List;

public enum TurboIDExperimentType {
	MIX, EMD(true), LBR(true), SUN1(true), MAN1(true), TURBOID_ONLY;
	private final boolean isBait;

	private TurboIDExperimentType() {
		this(false);
	}

	private TurboIDExperimentType(boolean isBait) {
		this.isBait = isBait;
	}

	public boolean isBait() {
		return this.isBait;
	}

	public static List<TurboIDExperimentType> getBaits() {
		final List<TurboIDExperimentType> ret = new ArrayList<TurboIDExperimentType>();
		for (final TurboIDExperimentType bait : values()) {
			if (bait.isBait) {
				ret.add(bait);
			}
		}
		return ret;
	}

	public static List<TurboIDExperimentType> getBaitsAndTBID() {
		final List<TurboIDExperimentType> baits = getBaits();
		baits.add(TurboIDExperimentType.TURBOID_ONLY);
		return baits;
	}
}
