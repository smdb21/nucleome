package edu.scripps.yates.nucleome;

import java.io.IOException;
import java.util.Collection;
import java.util.List;

import edu.scripps.yates.nucleome.model.CellCompartment;
import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.nucleome.model.Wash;

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

	public abstract double getScore(Collection<String> proteinAccessions, CellType celltype,
			CellCompartment cellCompartmentToStudy) throws IOException;

	public abstract double getScore(Collection<String> proteinAccessions, CellType celltype, Wash wash,
			CellCompartment cellCompartmentToStudy) throws IOException;

	public abstract double getScore(Collection<String> proteinAccessions, CellType celltype, Wash wash,
			CellCompartment cellCompartmentToStudy, Collection<CellCompartment> fractionTypes) throws IOException;

	public abstract double getScore(Collection<String> proteinAccessions, CellType celltype, Collection<Wash> washes,
			CellCompartment cellCompartmentToStudy, Collection<CellCompartment> fractionTypes) throws IOException;

	public abstract String getName();

	public abstract double getSumNSAFs(List<String> proteinAccessions, CellType cellType, Wash washPreWash,
			CellCompartment fraction) throws IOException;
}
