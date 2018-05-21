package edu.scripps.yates.nucleome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Date;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.nucleome.model.Wash;
import edu.scripps.yates.utilities.venndata.VennData;

public class TripleComparisonReport extends PairComparisonReport {

	private final Map<String, Double> scoreMap3;
	private final CellType cellType3;
	private final Wash wash3;

	public TripleComparisonReport(_4DNucleomeAnalyzer nucleomeAnalyzer, CellType cellType1, Wash wash1,
			Map<String, Double> scoreMap1, CellType cellType2, Wash wash2, Map<String, Double> scoreMap2,
			CellType cellType3, Wash wash3, Map<String, Double> scoreMap3) {
		super(nucleomeAnalyzer, cellType1, wash1, scoreMap1, cellType2, wash2, scoreMap2);

		this.scoreMap3 = scoreMap3;
		this.cellType3 = cellType3;
		this.wash3 = wash3;
	}

	@Override
	public VennData getVennData() throws IOException {

		if (venn == null) {
			Set<String> enriched1 = getEnriched(scoreMap1, cellType1, wash1);
			Set<String> enriched2 = getEnriched(scoreMap2, cellType2, wash2);
			Set<String> enriched3 = getEnriched(scoreMap3, cellType3, wash3);
			venn = new VennData(cellType1.name() + " vs " + cellType2.name() + " vs " + cellType3.name(),
					cellType1.name(), enriched1, cellType2.name(), enriched2, cellType3.name(), enriched3);
		}
		return venn;
	}

	@Override
	public void printToFile(File file) throws IOException {

		// pritn to the file
		FileWriter fw = null;
		try {
			fw = new FileWriter(file);
			fw.write("Report made on " + new Date() + "\n");
			fw.write("We consider enriched proteins when having an enrichment score > "
					+ Constants.ENRICHMENT_SCORE_THRESHOLD + "\n\n");
			fw.write(getEnriched(scoreMap1, cellType1, wash1).size() + " proteins enriched in "
					+ Constants.cellCompartmentToStudy + " in " + cellType1 + "\n");
			fw.write(getEnriched(scoreMap2, cellType2, wash2).size() + " proteins enriched in "
					+ Constants.cellCompartmentToStudy + " in " + cellType2 + "\n");
			fw.write(getEnriched(scoreMap3, cellType3, wash3).size() + " proteins enriched in "
					+ Constants.cellCompartmentToStudy + " in " + cellType3 + "\n");

			fw.write(getVennData().getUniqueTo1().size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in " + cellType1 + " and not in " + cellType2 + " and not in " + cellType3 + "\n");
			fw.write(getVennData().getUniqueTo2().size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in " + cellType2 + " and not in " + cellType1 + " and not in " + cellType3 + "\n");
			fw.write(getVennData().getUniqueTo3().size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in " + cellType3 + " and not in " + cellType1 + " and not in " + cellType2 + "\n");
			fw.write(getVennData().getIntersection123().size() + " proteins enriched in "
					+ Constants.cellCompartmentToStudy + " in all " + cellType1 + " and " + cellType2 + " and "
					+ cellType3 + "\n\n");

			fw.write("\nCopy and paste this URL to see the graphical representation of the comparison:\n"
					+ getVennData().getImageURL().toString() + "\n\n");
			fw.write("\n" + getVennData().getUniqueTo1().size() + " proteins enriched in "
					+ Constants.cellCompartmentToStudy + " in " + cellType1 + " and not in " + cellType2
					+ " and not in " + cellType3 + ":\n");
			fw.write("#\tACC\tEnrichment_Score\tgene name(s)\tprotein description(s)\n");
			writeListOfPairs(fw, getVennData().getUniqueTo1(), scoreMap1);

			fw.write("\n" + getVennData().getUniqueTo2().size() + " proteins enriched in "
					+ Constants.cellCompartmentToStudy + " in " + cellType2 + " and not in " + cellType1
					+ " and not in " + cellType3 + ":\n");
			fw.write("#\tACC\tEnrichment_Score\tprotein description\n");
			writeListOfPairs(fw, getVennData().getUniqueTo2(), scoreMap2);

			fw.write("\n" + getVennData().getUniqueTo3().size() + " proteins enriched in "
					+ Constants.cellCompartmentToStudy + " in " + cellType3 + " and not in " + cellType1
					+ " and not in " + cellType2 + ":\n");
			fw.write("#\tACC\tEnrichment_Score\tprotein description\n");
			writeListOfPairs(fw, getVennData().getUniqueTo3(), scoreMap3);

			fw.write("\nProteins enriched in " + Constants.cellCompartmentToStudy + " in all " + cellType1 + " and in "
					+ cellType2 + " and in " + cellType3 + ":\n");
			fw.write("#\tACC\tDelta_Enrichment_Score (" + cellType1 + "|" + cellType2 + "|" + cellType3
					+ ")\tprotein description\n");
			writeListOfPairs3(fw, getVennData().getIntersection123(), scoreMap1, scoreMap2, scoreMap3);

		} finally {
			fw.close();
		}

	}

	protected void writeListOfPairs3(FileWriter fw, Collection<Object> proteinAccs, Map<String, Double> scoreMap1,
			Map<String, Double> scoreMap2, Map<String, Double> scoreMap3) throws IOException {
		int i = 1;
		for (Object obj : proteinAccs) {
			String proteinAcc = obj.toString();

			Double score1 = scoreMap1.get(proteinAcc);
			Double score2 = scoreMap2.get(proteinAcc);
			Double score3 = scoreMap3.get(proteinAcc);
			String proteinNameString = nucleomeAnalyzer.getProteinNameString(proteinAcc, null, null);
			String geneNameString = nucleomeAnalyzer.getGeneNameString(proteinAcc, null, null);
			fw.write(i++ + "\t" + proteinAcc + "\t" + score1 + "|" + score2 + "|" + score3 + "\t" + geneNameString
					+ "\t" + proteinNameString + "\n");
		}
	}

}
