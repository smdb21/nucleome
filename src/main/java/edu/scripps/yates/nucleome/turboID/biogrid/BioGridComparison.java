package edu.scripps.yates.nucleome.turboID.biogrid;

import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.excel.ExcelReader;
import gnu.trove.map.hash.THashMap;

/**
 * This class will compare the known interactors in BioGrid for the baits
 * EMD-LBR-SUN1-MAN1 with the interactors discovered in our experiment
 * 
 * @author salvador
 *
 */
public class BioGridComparison {
	private static final Logger log = Logger.getLogger(BioGridComparison.class);
	private final static File biogridFolder = new File(
			"Z:\\share\\Salva\\data\\4D_Nucleome\\TurboID\\input\\optimized_params\\BioGrid");
	private final static File experimentalDataFile = new File(
			"Z:\\share\\Salva\\data\\4D_Nucleome\\TurboID\\input\\optimized_params\\Nu_anova_distribInt_ANOVA_TTESTs_DATA.xlsx");
	private final double minExperimentalScoreThreshold = 0.0;
	private final double maxExperimentalScoreThreshold = 0.40;
	private final double experimentalScoreThresholdIncrement = 0.01;

	private final String[] baits = { "EMD", "LBR", "MAN1", "SUN1" };

	public static void main(String[] args) {
		final BioGridComparison biogridComparator = new BioGridComparison();
		try {
			biogridComparator.run();
		} catch (final IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		System.out.println("Everything correct");
		System.exit(0);
	}

	private void run() throws IOException {
		final Map<String, BioGridInteractors> biogridInteractorsByBait = readBioGridInteractorsByBait();
		final Map<String, BioGridInteractors> preysByBait = readExperimentalPreysByBait();
		final FileWriter fw = new FileWriter(
				biogridFolder.getAbsolutePath() + File.separator + "BioGRID_comparison.txt");
		fw.write("bait" + "\t" + "numInCommon" + "\t" + "score threshold" + "\t" + "BIOGRID" + "\t"
				+ "experimental preys" + "\n");

		for (final String bait : baits) {
			double threshold = minExperimentalScoreThreshold;
			while (threshold <= maxExperimentalScoreThreshold) {

				final BioGridInteractors bioGridInteractors = biogridInteractorsByBait.get(bait);
				final BioGridInteractors preysObj = preysByBait.get(bait);
				final List<String> biogrid = bioGridInteractors.getInteractorsAsList();
				final Set<String> preys = preysObj.getInteractorByMinScore(threshold);
				final int numInCommon = getNumInCommon(biogrid, preys);
				fw.write(bait + "\t" + numInCommon + "\t" + threshold + "\t" + biogrid.size() + "\t" + preys.size()
						+ "\n");
				threshold += experimentalScoreThresholdIncrement;
			}

		}
		fw.close();
	}

	private int getNumInCommon(Collection<String> biogrid, Collection<String> preys) {
		int count = 0;
		for (final String bio : biogrid) {
			if (preys.contains(bio)) {
				count++;
			}
		}
		return count;
	}

	private Map<String, BioGridInteractors> readExperimentalPreysByBait() throws IOException {
		final Map<String, BioGridInteractors> ret = new THashMap<String, BioGridInteractors>();
		final ExcelReader reader = new ExcelReader(experimentalDataFile, 0, 0);
		int numRow = 1;
		final int emdIndex = reader.getColumnIndex(0, "distNormInt EMD");
		ret.put("EMD", new BioGridInteractors("EMD"));
		final int lbrIndex = reader.getColumnIndex(0, "distNormInt LBR");
		ret.put("LBR", new BioGridInteractors("LBR"));
		final int sun1Index = reader.getColumnIndex(0, "distNormInt SUN1");
		ret.put("SUN1", new BioGridInteractors("SUN1"));
		final int man1Index = reader.getColumnIndex(0, "distNormInt MAN1");
		ret.put("MAN1", new BioGridInteractors("MAN1"));
		final int geneIndex = reader.getColumnIndex(0, "gene");
		while (true) {
			final String gene = reader.getStringValue(0, numRow, geneIndex);
			if (gene == null) {
				break;
			}
			final String scoreEMD = reader.getNumberValue(0, numRow, emdIndex);
			if (!scoreEMD.equals("NaN")) {
				ret.get("EMD").addInteractor(gene, Double.valueOf(scoreEMD));
			}
			final String scoreLBR = reader.getNumberValue(0, numRow, lbrIndex);
			if (!scoreLBR.equals("NaN")) {
				ret.get("LBR").addInteractor(gene, Double.valueOf(scoreLBR));
			}
			final String scoreSUN1 = reader.getNumberValue(0, numRow, sun1Index);
			if (!scoreSUN1.equals("NaN")) {
				ret.get("SUN1").addInteractor(gene, Double.valueOf(scoreSUN1));
			}
			final String scoreMAN1 = reader.getNumberValue(0, numRow, man1Index);
			if (!scoreMAN1.equals("NaN")) {
				ret.get("MAN1").addInteractor(gene, Double.valueOf(scoreMAN1));
			}
			numRow++;
		}
		log.info(ret.get("EMD").getInteractors().size() + " interactors for EMD");
		log.info(ret.get("LBR").getInteractors().size() + " interactors for LBR");
		log.info(ret.get("SUN1").getInteractors().size() + " interactors for SUN1");
		log.info(ret.get("MAN1").getInteractors().size() + " interactors for MAN1");
		return ret;
	}

	private Map<String, BioGridInteractors> readBioGridInteractorsByBait() throws IOException {
		final File[] biogridFiles = biogridFolder.listFiles(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String name) {
				if (name.endsWith("tab2.txt")) {
					return true;
				}
				return false;
			}
		});
		final Map<String, BioGridInteractors> ret = new THashMap<String, BioGridInteractors>();
		for (final File biogridFile : biogridFiles) {
			final String fileName = FilenameUtils.getBaseName(biogridFile.getAbsolutePath());
			final String bait = fileName.substring(0, fileName.indexOf("_"));
			final BioGridInteractors interactors = readBioGridFile(biogridFile, bait);
			ret.put(bait, interactors);
			log.info(interactors.getInteractors().size() + " for " + bait + " in BIOGRID");
		}

		return ret;
	}

	private BioGridInteractors readBioGridFile(File biogridFile, String bait) throws IOException {
		final BioGridInteractors ret = new BioGridInteractors("BIOGRID");
		final List<String> lines = readAllLinesFromZip(biogridFile);
		final Map<String, Integer> indexesByHEaders = getIndexesByHeaders(lines.get(0));
		for (int numRow = 1; numRow < lines.size(); numRow++) {
			final String line = lines.get(numRow);
			final int numCol1 = indexesByHEaders.get("Official Symbol Interactor A");
			final int numCol2 = indexesByHEaders.get("Official Symbol Interactor B");
			final String[] split = line.split("\t");
			final String interactor1 = split[numCol1];
			final String interactor2 = split[numCol2];
			if (!interactor1.equals(bait)) {
				ret.addInteractor(interactor1, null);
			} else {
				ret.addInteractor(interactor2, null);
			}
		}
		return ret;
	}

	private List<String> readAllLinesFromZip(File file) throws IOException {
		return Files.readAllLines(file.toPath());
	}

	private Map<String, Integer> getIndexesByHeaders(String wholeHeader) {
		final Map<String, Integer> ret = new THashMap<String, Integer>();
		final String[] split = wholeHeader.split("\t");
		int col = 0;
		for (final String header : split) {
			ret.put(header, col);
			col++;
		}
		return ret;
	}
}
