package edu.scripps.yates.nucleome;

import java.io.IOException;

import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.utilities.grouping.ProteinGroup;

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

	public abstract double getScore(ProteinGroup proteinGroup, CellType celltype) throws IOException;
}
