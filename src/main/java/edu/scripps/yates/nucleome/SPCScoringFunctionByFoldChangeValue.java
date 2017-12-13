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
 * This is the second scoring function proposed by LArry. It is based on the SPC
 * difference but not the ratio.
 * 
 * @author Salva
 *
 */
public class SPCScoringFunctionByFoldChangeValue extends ScoringFunction {
	private final static Logger log = Logger.getLogger(SPCScoringFunctionByFoldChangeValue.class);
	private final _4DNucleomeAnalyzer analyzer;

	public SPCScoringFunctionByFoldChangeValue(_4DNucleomeAnalyzer analyzer) {
		this.analyzer = analyzer;
	}

	@Override
	public double getScore(Collection<String> proteinAccessions, CellType celltype) throws IOException {
		int numExperiments = 0;
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

		// log.info(experimentList.size() + " experiments");
		for (Experiment experiment : experimentList) {
			double experimentScore = 0.0;
			int numComparisons = 0;
			boolean ratioInThisExperiment = false;

			for (CellCompartment denominatorCellCompartment : CellCompartment.values()) {
				if (denominatorCellCompartment == Constants.cellCompartmentToStudy) {
					continue;
				}
				double spcToStudy = experiment.getAvgSpectralCount(proteinAccessions, Constants.cellCompartmentToStudy,
						true);
				double spcDenominator = experiment.getAvgSpectralCount(proteinAccessions, denominatorCellCompartment,
						true);
				if (spcToStudy > 0 || spcDenominator > 0) {
					ratioInThisExperiment = true;
					numComparisons++;
				} else if (spcToStudy == 0 && spcDenominator == 0) {
					continue;
				}
				if (spcToStudy >= spcDenominator) {
					experimentScore += 1;
					if (spcToStudy >= 3 * spcDenominator && spcToStudy > 3) {
						experimentScore += 1;
					}
				} else {
					if (Constants.includeNegativeScoring) {
						experimentScore -= 1;
						if (3 * spcToStudy <= spcDenominator && spcDenominator > 3) {
							experimentScore -= 1;
						}
					}
				}
			}
			// normalize by the number of comparisons
			if (numComparisons > 0) {
				experimentScore = experimentScore / numComparisons;
			}
			score += experimentScore;
			if (ratioInThisExperiment) {
				numExperiments++;
			}

		}
		// if we are looking in a particular cellType, divide
		score = score / numExperiments;
		return score;
	}
}
