package edu.scripps.yates.nucleome;

import java.io.IOException;
import java.util.Collection;
import java.util.List;

import org.apache.log4j.Logger;

import edu.scripps.yates.nucleome.model.CellCompartment;
import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.nucleome.model.Experiment;
import edu.scripps.yates.nucleome.model.Wash;

/**
 * 
 * 
 * @author Salva
 *
 */
public class ScoringFunctionByNE_NSAF_Ratios extends ScoringFunction {
	private final static Logger log = Logger.getLogger(ScoringFunctionByNE_NSAF_Ratios.class);
	private final _4DNucleomeAnalyzer analyzer;

	public ScoringFunctionByNE_NSAF_Ratios(_4DNucleomeAnalyzer analyzer) {
		this.analyzer = analyzer;
	}

	@Override
	public double getScore(Collection<String> proteinAccessions, CellType celltype,
			CellCompartment cellCompartmentToStudy) throws IOException {
		return getScore(proteinAccessions, celltype, null, cellCompartmentToStudy);
	}

	@Override
	public double getScore(Collection<String> proteinAccessions, CellType celltype, Wash wash,
			CellCompartment cellCompartmentToStudy) throws IOException {
		double numerator = 0.0;
		double denominator = 0.0;
		int totalSPC = 0;

		final List<Experiment> experimentList = analyzer.getExperiments(celltype, wash);
		if (proteinAccessions.contains("Q8K3Z9")) {
			log.info(proteinAccessions);
		}
		// log.info(experimentList.size() + " experiments");
		for (final Experiment experiment : experimentList) {

			final double numNSAF = experiment.getSumNSAF(proteinAccessions, cellCompartmentToStudy, true);
			numerator += numNSAF;
			denominator += numNSAF;
			final int totalSPCExperiment = experiment.getSumSPC(proteinAccessions, cellCompartmentToStudy, true);
			totalSPC += totalSPCExperiment;
			for (final CellCompartment denominatorCellCompartment : CellCompartment.values()) {
				if (denominatorCellCompartment == cellCompartmentToStudy) {
					continue;
				}
				final double sumNSAF = experiment.getSumNSAF(proteinAccessions, denominatorCellCompartment, true);
				denominator += sumNSAF;
			}

		}

		if (denominator > 0.0) {
			final double ratio = numerator / denominator;
			return ratio;
		} else {
			// it depends on the number in the numerator
			if (totalSPC > Constants.MIN_TOTAL_SPC) {
				return Double.POSITIVE_INFINITY;
			} else {
				return Double.NaN;
			}
		}

	}

	@Override
	public double getScore(Collection<String> proteinAccessions, CellType celltype, Wash wash,
			CellCompartment cellCompartmentToStudy, Collection<CellCompartment> fractionTypes) throws IOException {
		double numerator = 0.0;
		double denominator = 0.0;

		final List<Experiment> experimentList = analyzer.getExperiments(celltype, wash);
		// log.info(experimentList.size() + " experiments");
		for (final Experiment experiment : experimentList) {

			final double numNSAF = experiment.getSumNSAF(proteinAccessions, cellCompartmentToStudy, true);
			numerator += numNSAF;
			denominator += numNSAF;
			for (final CellCompartment denominatorCellCompartment : CellCompartment.values()) {
				if (denominatorCellCompartment == cellCompartmentToStudy
						// only the ones in the collection of fractionTypes are
						// considered
						|| (fractionTypes != null && !fractionTypes.contains(denominatorCellCompartment))) {
					continue;
				}
				final double sumNSAF = experiment.getSumNSAF(proteinAccessions, denominatorCellCompartment, true);
				denominator += sumNSAF;
			}

		}

		// if (denominator > 0.0) {
		final double ratio = numerator / denominator;
		if (Double.isInfinite(ratio)) {
			log.info(numerator + "\t" + denominator);
		}
		return ratio;
		// } else {
		// // it depends on the number in the numerator
		// if (totalSPC > Constants.MIN_TOTAL_SPC) {
		// return Double.POSITIVE_INFINITY;
		// } else {
		// return Double.NaN;
		// }
		// }

	}

	@Override
	public String getName() {
		return "NSAFratios";
	}

	@Override
	public double getScore(Collection<String> proteinAccessions, CellType celltype, Collection<Wash> washes,
			CellCompartment cellCompartmentToStudy, Collection<CellCompartment> fractionTypes) throws IOException {
		double numerator = 0.0;
		double denominator = 0.0;
		int totalSPC = 0;

		final List<Experiment> experimentList = analyzer.getExperiments(celltype, washes);
		// log.info(experimentList.size() + " experiments");
		for (final Experiment experiment : experimentList) {

			final double numNSAF = experiment.getSumNSAF(proteinAccessions, cellCompartmentToStudy, true);
			numerator += numNSAF;
			denominator += numNSAF;
			final int totalSPCExperiment = experiment.getSumSPC(proteinAccessions, cellCompartmentToStudy, true);
			totalSPC += totalSPCExperiment;
			for (final CellCompartment denominatorCellCompartment : CellCompartment.values()) {
				if (denominatorCellCompartment == cellCompartmentToStudy
						// only the ones in the collection of fractionTypes are
						// considered
						|| (fractionTypes != null && !fractionTypes.contains(denominatorCellCompartment))) {
					continue;
				}
				final double sumNSAF = experiment.getSumNSAF(proteinAccessions, denominatorCellCompartment, true);
				denominator += sumNSAF;
			}

		}

		if (denominator > 0.0) {
			final double ratio = numerator / denominator;
			return ratio;
		} else {
			// it depends on the number in the numerator
			if (totalSPC > Constants.MIN_TOTAL_SPC) {
				return Double.POSITIVE_INFINITY;
			} else {
				return Double.NaN;
			}
		}

	}

	@Override
	public double getSumNSAFs(List<String> proteinAccessions, CellType cellType, Wash washPreWash,
			CellCompartment fraction) throws IOException {
		double nsaf = 0.0;
		final List<Experiment> experimentList = analyzer.getExperiments(cellType, washPreWash);
		for (final Experiment experiment : experimentList) {
			final double numNSAF = experiment.getSumNSAF(proteinAccessions, fraction, true);
			nsaf += numNSAF;
		}
		return nsaf;
	}
}
