package edu.scripps.yates.nucleome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Date;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.utilities.venndata.VennData;

public class PairComparisonReport {
	protected final _4DNucleomeAnalyzer nucleomeAnalyzer;
	protected final Map<String, Double> scoreMap1;
	protected final Map<String, Double> scoreMap2;
	protected final CellType cellType1;
	protected final CellType cellType2;

	protected VennData venn;

	public PairComparisonReport(_4DNucleomeAnalyzer nucleomeAnalyzer, CellType cellType1, Map<String, Double> scoreMap1,
			CellType cellType2, Map<String, Double> scoreMap2) {
		this.nucleomeAnalyzer = nucleomeAnalyzer;
		this.scoreMap1 = scoreMap1;
		this.scoreMap2 = scoreMap2;
		this.cellType1 = cellType1;
		this.cellType2 = cellType2;
	}

	public VennData getVennData() throws IOException {

		if (venn == null) {
			Set<String> enriched1 = getEnriched(scoreMap1);
			Set<String> enriched2 = getEnriched(scoreMap2);
			venn = new VennData(cellType1.name() + " vs " + cellType2.name(), cellType1.name(), enriched1,
					cellType2.name(), enriched2, null, null);
		}
		return venn;
	}

	protected Set<String> getEnriched(Map<String, Double> scoreMap) throws IOException {
		Set<String> ret = new HashSet<String>();
		for (String proteinAcc : scoreMap.keySet()) {
			Double enrichmentScore = scoreMap.get(proteinAcc);
			if (!this.nucleomeAnalyzer.isValid(proteinAcc)) {
				continue;
			}
			if (!Double.isNaN(enrichmentScore) && enrichmentScore >= Constants.ENRICHMENT_SCORE_THRESHOLD) {
				ret.add(proteinAcc);
			}
		}
		return ret;
	}

	public void printToFile(File file) throws IOException {

		// pritn to the file
		FileWriter fw = null;
		try {
			fw = new FileWriter(file);
			fw.write("Report made on " + new Date() + "\n");
			fw.write("We consider enriched proteins when having an enrichment score > "
					+ Constants.ENRICHMENT_SCORE_THRESHOLD + "\n\n");
			fw.write(getEnriched(scoreMap1).size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in " + cellType1 + "\n");
			fw.write(getEnriched(scoreMap2).size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in " + cellType2 + "\n");
			VennData vennData = getVennData();
			fw.write(vennData.getUniqueTo1().size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in " + cellType1 + " and not in " + cellType2 + "\n");
			fw.write(vennData.getUniqueTo2().size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in " + cellType2 + " and not in " + cellType1 + "\n");
			fw.write(vennData.getIntersection12().size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in both " + cellType1 + " and " + cellType2 + "\n\n");

			fw.write("\nCopy and paste this URL to see the graphical representation of the comparison:\n"
					+ vennData.getImageURL().toString() + "\n\n");

			fw.write("\n" + vennData.getUniqueTo1().size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in " + cellType2 + " and not in " + cellType1 + ":\n");
			fw.write("#\tACC(s)\tEnrichment_Score\tgene name(s)\tprotein description(s)\n");
			writeListOfPairs(fw, vennData.getUniqueTo1(), scoreMap1);

			fw.write("\n" + vennData.getUniqueTo2().size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in " + cellType2 + " and not in " + cellType1 + ":\n");
			fw.write("#\tACC\tEnrichment_Score\tprotein description\n");
			writeListOfPairs(fw, vennData.getUniqueTo2(), scoreMap2);

			fw.write("\n" + vennData.getIntersection12().size() + " proteins enriched in "
					+ Constants.cellCompartmentToStudy + " in both " + cellType1 + " and in " + cellType2 + ":\n");
			fw.write("#\tACC\tDelta_Enrichment_Score (" + cellType1 + "|" + cellType2 + ")\tprotein description\n");
			writeListOfPairs2(fw, vennData.getIntersection12(), scoreMap1, scoreMap2);

		} finally {
			fw.close();
		}

	}

	protected void writeListOfPairs(FileWriter fw, Collection<Object> proteinAccs, Map<String, Double> scoreMap)
			throws IOException {
		int i = 1;

		for (Object obj : proteinAccs) {
			String proteinAcc = obj.toString();
			Double score = scoreMap.get(proteinAcc);

			String proteinNameString = nucleomeAnalyzer.getProteinNameString(proteinAcc);
			String geneNameString = nucleomeAnalyzer.getGeneNameString(proteinAcc);
			fw.write(i++ + "\t" + proteinAcc + "\t" + score + "\t" + geneNameString + "\t" + proteinNameString + "\n");
		}
	}

	protected double getSumScores(String secondElement) {
		double ret = 0.0;
		final String[] split = secondElement.split("\\|");
		for (String string : split) {
			ret += Double.valueOf(string.trim());
		}
		return ret;
	}

	protected void writeListOfPairs2(FileWriter fw, Collection<Object> proteinAccs, Map<String, Double> scoreMap1,
			Map<String, Double> scoreMap2) throws IOException {
		int i = 1;
		for (Object obj : proteinAccs) {
			String proteinAcc = obj.toString();

			Double score1 = scoreMap1.get(proteinAcc);
			Double score2 = scoreMap2.get(proteinAcc);
			String proteinNameString = nucleomeAnalyzer.getProteinNameString(proteinAcc);
			String geneNameString = nucleomeAnalyzer.getGeneNameString(proteinAcc);
			fw.write(i++ + "\t" + proteinAcc + "\t" + score1 + "|" + score2 + "\t" + geneNameString + "\t"
					+ proteinNameString + "\n");
		}
	}
}
