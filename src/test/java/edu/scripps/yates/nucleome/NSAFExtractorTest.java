package edu.scripps.yates.nucleome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FilenameUtils;
import org.junit.Test;

import edu.scripps.yates.annotations.uniprot.UniprotProteinRetriever;
import edu.scripps.yates.nucleome.model.CellCompartment;
import edu.scripps.yates.nucleome.model.Experiment;
import edu.scripps.yates.nucleome.model.Fractionation;
import edu.scripps.yates.nucleome.model.Replicate;

public class NSAFExtractorTest {
	@Test
	public void extractNSAFValues() {
		String[] genes = { "Lmna", "Lmnb1", "Lmnb2", "Tmpo", "Emd", "Lemd3", "Lemd2", "Tor1aip1", "Lbr", "Sun1", "Sun2",
				"Syne1", "Syne2", "Syne3", "Prr14", "Pom121", "Ndc1", "Nup160", "Nup133", "Nup107", "Nup85", "Nup43",
				"Nup37", "Nup155", "Nup188", "Nup93", "Nup35" };
		_4DNucleomeAnalyzer analyzer;
		FileWriter fw = null;
		try {

			String pass = null;
			analyzer = new _4DNucleomeAnalyzer(pass);

			////////////////////////////////////////////////////////////
			// PARAMETERS
			Constants.GO_FILTER = false;
			Constants.cellCompartmentToStudy = CellCompartment.NE;
			Constants.TESTING = false;
			Constants.DATASET_PATHS_FILE = "z:\\share\\Salva\\data\\4D_Nucleome\\datasets_paths.txt";
			Constants.MIN_TOTAL_SPC = 5;
			Constants.geneFilter = genes;
			Constants.printScoreDistributions = false;
			Constants.compareScores = false;
			UniprotProteinRetriever.enableCache = false;

			////////////////////////////////////////////////////////////
			// load the data
			analyzer.run();
			fw = new FileWriter(new File(FilenameUtils.getFullPath(Constants.DATASET_PATHS_FILE) + File.separator
					+ "NSAF_NE_over_differentiation.txt"));
			// header
			StringBuilder sb2 = new StringBuilder("Gene\t");
			List<Experiment> experiments = analyzer.getAllExperiments();

			for (Experiment experiment : experiments) {
				List<Replicate> replicates = experiment.getReplicates();
				for (Replicate replicate : replicates) {
					Fractionation fractionation = replicate.getFractionation(CellCompartment.NE);
					sb2.append(fractionation.getName() + "\t");
				}
			}
			fw.write(sb2.toString() + "\n");
			// iterate over genes of interest
			for (String gene : genes) {

				StringBuilder sb = new StringBuilder(gene + "\t");

				for (Experiment experimentU : experiments) {
					List<Replicate> replicates = experimentU.getReplicates();
					for (Replicate replicate : replicates) {
						Fractionation fractionation = replicate.getFractionation(CellCompartment.NE);
						Set<String> proteinAccessions = fractionation.getProteins().stream()
								.filter(p -> p.getGenes().stream().anyMatch(g -> g.getGeneID().equalsIgnoreCase(gene)))
								.map(p -> p.getAccession()).collect(Collectors.toSet());
						if (proteinAccessions != null && !proteinAccessions.isEmpty()) {
							double sumNSAF = fractionation.getAverageNSAF(proteinAccessions, false) * 1000000;
							sb.append(sumNSAF);
						}
						sb.append("\t");
					}
				}
				fw.write(sb.toString() + "\n");
			}
			System.out.println("DONE");

		} catch (RuntimeException e) {
			e.printStackTrace();
			System.err.println("ERROR: " + e.getMessage());
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("ERROR: " + e.getMessage());
		} finally {
			if (fw != null) {
				try {
					fw.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

	}
}
