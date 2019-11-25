package edu.scripps.yates.nucleome.turboID;

import java.util.ArrayList;
import java.util.List;

public enum Replicate {
	A(TurboIDFraction.NU), B(TurboIDFraction.NU), Brerun(TurboIDFraction.NU), Arerun(TurboIDFraction.NU),
	CyA(TurboIDFraction.CY), CyB(TurboIDFraction.CY);
	private final TurboIDFraction fraction;

	private Replicate(TurboIDFraction fraction) {
		this.fraction = fraction;
	}

	public static Replicate[] values(TurboIDFraction fraction) {
		final List<Replicate> ret = new ArrayList<Replicate>();
		for (final Replicate replicate : values()) {
			if (replicate.getFraction() == fraction) {
				ret.add(replicate);
			}
		}
		return ret.toArray(new Replicate[0]);
	}

	public TurboIDFraction getFraction() {
		return fraction;
	}
}
