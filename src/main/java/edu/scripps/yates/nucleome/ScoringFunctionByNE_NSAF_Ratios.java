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
	public double getScore(Collection<String> proteinAccessions, CellType celltype) throws IOException {
		return getScore(proteinAccessions, celltype, null);
	}

	@Override
	public double getScore(Collection<String> proteinAccessions, CellType celltype, Wash wash) throws IOException {
		double numerator = 0.0;
		double denominator = 0.0;
		int totalSPC = 0;

		List<Experiment> experimentList = analyzer.getExperiments(celltype, wash);
		if (proteinAccessions.contains("Q8K3Z9")) {
			log.info(proteinAccessions);
		}
		// log.info(experimentList.size() + " experiments");
		for (Experiment experiment : experimentList) {

			double numNSAF = experiment.getSumNSAF(proteinAccessions, Constants.cellCompartmentToStudy, true);
			numerator += numNSAF;
			denominator += numNSAF;
			int totalSPCExperiment = experiment.getSumSPC(proteinAccessions, Constants.cellCompartmentToStudy, true);
			totalSPC += totalSPCExperiment;
			for (CellCompartment denominatorCellCompartment : CellCompartment.values()) {
				if (denominatorCellCompartment == Constants.cellCompartmentToStudy) {
					continue;
				}
				double sumNSAF = experiment.getSumNSAF(proteinAccessions, denominatorCellCompartment, true);
				denominator += sumNSAF;
			}

		}

		if (denominator > 0.0) {
			double ratio = numerator / denominator;
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
	public String getName() {
		return "NSAFratios";
	}
}
