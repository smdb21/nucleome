package edu.scripps.yates.nucleome;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.log4j.Logger;

import edu.scripps.yates.nucleome.model.CellCompartment;
import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.nucleome.model.Experiment;

/**
 * Calculates the % of SPC that are in NE over N and C
 * 
 * @author Salva
 *
 */
public class ScoringFunctionByNE_NSAF_Percentage extends ScoringFunction {
	private final static Logger log = Logger.getLogger(ScoringFunctionByNE_NSAF_Percentage.class);
	private final _4DNucleomeAnalyzer analyzer;

	public ScoringFunctionByNE_NSAF_Percentage(_4DNucleomeAnalyzer analyzer) {
		this.analyzer = analyzer;
	}

	@Override
	public double getScore(Collection<String> proteinAccessions, CellType celltype) throws IOException {
		double score = 0.0;
		List<Experiment> experimentList = new ArrayList<Experiment>();
		if (celltype == null || celltype == CellType.A) {
			experimentList.addAll(analyzer.getExperimentsA());
		}
		if (celltype == null || celltype == CellType.M) {
			experimentList.addAll(analyzer.getExperimentsM());
		}
		if (celltype == null || celltype == CellType.U) {
			experimentList.addAll(analyzer.getExperimentsU());
		}
		int numExperimentsInWhichIsDetected = 0;
		// log.info(experimentList.size() + " experiments");
		for (Experiment experiment : experimentList) {
			boolean presentInThisExperiment = false;
			int sumSPC = experiment.getSumSPC(proteinAccessions, Constants.cellCompartmentToStudy, true);
			if (sumSPC > 0) {
				presentInThisExperiment = true;
			}
			double avgNSAF = experiment.getAvgNSAF(proteinAccessions, Constants.cellCompartmentToStudy, true);

			for (CellCompartment denominatorCellCompartment : CellCompartment.values()) {
				if (denominatorCellCompartment == Constants.cellCompartmentToStudy) {
					continue;
				}
				int sumSPCOther = experiment.getSumSPC(proteinAccessions, denominatorCellCompartment, true);
				if (!presentInThisExperiment && sumSPCOther > 0) {
					presentInThisExperiment = true;
				}
				double avgNSAFOther = experiment.getAvgNSAF(proteinAccessions, denominatorCellCompartment, true);
				if (avgNSAF > avgNSAFOther) {
					score = score + 1;
				}
				if (avgNSAF > 3 * avgNSAFOther) {
					score = score + 1;
				}
			}
			if (presentInThisExperiment) {
				numExperimentsInWhichIsDetected++;
			}
		}
		if (numExperimentsInWhichIsDetected == 0.0) {
			return Double.NaN;
		}
		// normalize score
		score = score / numExperimentsInWhichIsDetected;
		return score;
	}
}
