package edu.scripps.yates.nucleome;

import java.io.IOException;
import java.util.Collection;

import edu.scripps.yates.nucleome.model.CellType;

public abstract class ScoringFunction {
	public static double getLog2Ratio(double ratio) {
		if (Double.isNaN(ratio)) {
			return Double.NaN;
		}
		if (ratio == Double.POSITIVE_INFINITY) {
			return Double.POSITIVE_INFINITY;
		}
		if (ratio == 0.0) {
			return Double.NEGATIVE_INFINITY;
		}

		return Math.log(ratio) / Math.log(2);
	}

	// public abstract double getScore(String proteinAcc, CellType celltype)
	// throws IOException;

	public abstract double getScore(Collection<String> proteinAccessions, CellType celltype) throws IOException;
}
