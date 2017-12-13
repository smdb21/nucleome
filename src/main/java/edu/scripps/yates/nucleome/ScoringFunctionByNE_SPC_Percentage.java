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
public class ScoringFunctionByNE_SPC_Percentage extends ScoringFunction {
	private final static Logger log = Logger.getLogger(ScoringFunctionByNE_SPC_Percentage.class);
	private final _4DNucleomeAnalyzer analyzer;

	public ScoringFunctionByNE_SPC_Percentage(_4DNucleomeAnalyzer analyzer) {
		this.analyzer = analyzer;
	}

	@Override
	public double getScore(Collection<String> proteinAccessions, CellType celltype) throws IOException {
		double score = 0.0;
		if (proteinAccessions.contains("E9QP46")) {
			log.info("aasdf ");
		}
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
		double spcToStudy = 0.0;
		double spcTotal = 0.0;
		// log.info(experimentList.size() + " experiments");
		for (Experiment experiment : experimentList) {
			int sumSPC = experiment.getSumSPC(proteinAccessions, Constants.cellCompartmentToStudy, true);
			spcToStudy += sumSPC;
			spcTotal += sumSPC;
			for (CellCompartment denominatorCellCompartment : CellCompartment.values()) {
				if (denominatorCellCompartment == Constants.cellCompartmentToStudy) {
					continue;
				}
				double spcOther = experiment.getSumSPC(proteinAccessions, denominatorCellCompartment, true);
				spcTotal += spcOther;
			}

		}
		if (spcTotal == 0.0) {
			return Double.NaN;
		}
		score = spcToStudy / spcTotal;
		return score;
	}
}
