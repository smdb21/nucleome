package edu.scripps.yates.nucleome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinRetriever;
import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.nucleome.filters.Filter;
import edu.scripps.yates.nucleome.filters.GOFilter;
import edu.scripps.yates.nucleome.model.CellCompartment;
import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.nucleome.model.ControlNE;
import edu.scripps.yates.nucleome.model.Experiment;
import edu.scripps.yates.utilities.dates.DatesUtil;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.PAnalyzer;
import edu.scripps.yates.utilities.grouping.ProteinEvidence;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.proteomicsmodel.Gene;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.ProteinAnnotation;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import edu.scripps.yates.utilities.util.Pair;

public class _4DNucleomeAnalyzer {
	private final static Logger log = Logger.getLogger(_4DNucleomeAnalyzer.class);

	public static void main(String[] args) {
		_4DNucleomeAnalyzer analyzer;
		try {
			analyzer = new _4DNucleomeAnalyzer();
			analyzer.run();

			System.out.println("DONE");
			System.exit(0);
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("ERROR: " + e.getMessage());
		}
		System.exit(-1);
	}

	private HashSet<String> proteinAccs;
	private final List<Experiment> experimentsU = new ArrayList<Experiment>();
	private final List<Experiment> experimentsA = new ArrayList<Experiment>();
	private final List<Experiment> experimentsM = new ArrayList<Experiment>();
	private ScoringFunction scoringFunction;
	private final String hostName = "jaina.scripps.edu";;
	private final String userName = "salvador";
	private final String pass = "Natjeija21";
	private final String remotefileName = "DTASelect-filter.txt";
	private final String datasetsPathsFile = "z:\\share\\Salva\\data\\4D_Nucleome\\datasets_paths.txt";
	private final String datasetsPhosphoPathsFile = "z:\\share\\Salva\\data\\4D_Nucleome\\datasets_paths_phospho.txt";
	private final File outputFolder = new File("z:\\share\\Salva\\data\\4D_Nucleome\\analysis");
	private final Map<CellType, List<Pair<String, Double>>> scoresByCellType = new HashMap<CellType, List<Pair<String, Double>>>();
	private List<ProteinGroup> proteinGroups;
	private GOFilter goFilter;

	public _4DNucleomeAnalyzer() throws IOException {

		////////////////////////////////////////////////////////////
		// PARAMETERS
		Constants.includeNegativeScoring = true;
		Constants.MIN_PEPTIDES_PER_PROTEIN = 1;
		Constants.MIN_PSM_PER_PROTEIN = 2;
		Constants.MIN_AVG_SPC = 3;
		Constants.GO_FILTER = true;
		Constants.cellCompartmentToStudy = CellCompartment.NE;
		Constants.TESTING = false;
		Constants.ENRICHMENT_SCORE_THRESHOLD = 3;
		Constants.DATASET_PATHS_FILE = datasetsPathsFile;
		UniprotProteinRetriever.enableCache = false;
		////////////////////////////////////////////////////////////

	}

	/**
	 * Returns a list of pairs (Protein ACC - Score), sorted by the score
	 *
	 * @return
	 * @throws IOException
	 */
	private List<Pair<String, Double>> calculateScoresFromGroups(CellType celltype) throws IOException {
		log.info("Calculating scores for " + proteinAccs.size() + " proteins in " + celltype);
		List<Pair<String, Double>> scores = new ArrayList<Pair<String, Double>>();
		for (ProteinGroup proteinGroup : getProteinGroups()) {
			if (proteinGroup.getKey().contains("D3Z008")) {
				System.out.println(proteinGroup);
			}
			double score = scoringFunction.getScore(proteinGroup, celltype);

			Pair<String, Double> pair = new Pair<String, Double>(proteinGroup.getKey(), score);
			scores.add(pair);
			// if (Constants.TESTING && !Double.isNaN(pair.getSecondElement()))
			// {
			// break;
			// }
		}
		log.info("Sorting scores...");
		// sort them by the score
		Collections.sort(scores, new Comparator<Pair<String, Double>>() {

			@Override
			public int compare(Pair<String, Double> o1, Pair<String, Double> o2) {
				final int compare = Double.compare(o2.getSecondElement(), o1.getSecondElement());
				return compare;
			}
		});
		log.info(scores.size() + " scores calculated");

		scoresByCellType.put(celltype, scores);
		return scores;
	}

	private List<ProteinGroup> getProteinGroups() throws IOException {
		if (proteinGroups == null) {
			PAnalyzer panalyzer = new PAnalyzer(false);
			Set<GroupableProtein> groupableProteins = getAllGroupableProteins();
			proteinGroups = panalyzer.run(groupableProteins);
		}
		return proteinGroups;
	}

	private Set<GroupableProtein> getAllGroupableProteins() throws IOException {
		Set<GroupableProtein> ret = new HashSet<GroupableProtein>();
		for (Experiment experiment : getExperimentsA()) {
			for (Protein protein : experiment.getProteins()) {
				ret.add(protein);
			}
		}
		for (Experiment experiment : getExperimentsM()) {
			for (Protein protein : experiment.getProteins()) {
				ret.add(protein);
			}
		}
		for (Experiment experiment : getExperimentsU()) {
			for (Protein protein : experiment.getProteins()) {
				ret.add(protein);
			}
		}
		return ret;
	}

	private void loadDatasets() throws IOException {
		log.info("Loading datasets");
		long t1 = System.currentTimeMillis();
		proteinAccs = new HashSet<String>();
		experimentsU.clear();
		experimentsA.clear();
		experimentsM.clear();

		DataPaths dataPaths = new DataPaths(Constants.DATASET_PATHS_FILE);
		// U (N, Ne, C)

		// FDR 1%

		String[] u2Files = dataPaths.getFiles("U2");
		String[] u3Files = dataPaths.getFiles("U3");
		String[] u4Files = dataPaths.getFiles("U41");
		String[] u42Files = dataPaths.getFiles("U42");
		String[] u5Files = dataPaths.getFiles("U5");
		String[] a2Files = dataPaths.getFiles("A2");
		String[] a4Files = dataPaths.getFiles("A41");
		String[] a42Files = dataPaths.getFiles("A42");
		String[] m1Files = dataPaths.getFiles("M11");
		String[] m12Files = dataPaths.getFiles("M12");
		String[] m3Files = dataPaths.getFiles("M31");
		String[] m32Files = dataPaths.getFiles("M32");

		Experiment experimentU2 = new Experiment("U2", CellType.U);
		experimentsU.add(experimentU2);
		experimentU2.addReplicate(1, CellType.U, CellCompartment.N, getRemoteFile(u2Files[0]));
		experimentU2.addReplicate(1, CellType.U, CellCompartment.NE, getRemoteFile(u2Files[1]));
		experimentU2.addReplicate(1, CellType.U, CellCompartment.C, getRemoteFile(u2Files[2]));
		proteinAccs.addAll(experimentU2.getProteinAccs());

		Experiment experimentU3 = new Experiment("U3", CellType.U);
		experimentsU.add(experimentU3);
		experimentU3.addReplicate(1, CellType.U, CellCompartment.N, getRemoteFile(u3Files[0]));
		experimentU3.addReplicate(1, CellType.U, CellCompartment.NE, getRemoteFile(u3Files[1]));
		experimentU3.addReplicate(1, CellType.U, CellCompartment.C, getRemoteFile(u3Files[2]));
		proteinAccs.addAll(experimentU3.getProteinAccs());
		Experiment experimentU4 = new Experiment("U4", CellType.U);
		experimentsU.add(experimentU4);
		experimentU4.addReplicate(1, CellType.U, CellCompartment.N, getRemoteFile(u4Files[0]));
		experimentU4.addReplicate(1, CellType.U, CellCompartment.NE, getRemoteFile(u4Files[1]));
		experimentU4.addReplicate(1, CellType.U, CellCompartment.C, getRemoteFile(u4Files[2]));
		experimentU4.addReplicate(2, CellType.U, CellCompartment.N, getRemoteFile(u42Files[0]));
		experimentU4.addReplicate(2, CellType.U, CellCompartment.NE, getRemoteFile(u42Files[1]));
		experimentU4.addReplicate(2, CellType.U, CellCompartment.C, getRemoteFile(u42Files[2]));
		proteinAccs.addAll(experimentU4.getProteinAccs());
		Experiment experimentU5 = new Experiment("U5", CellType.U);
		experimentsU.add(experimentU5);
		experimentU5.addReplicate(1, CellType.U, CellCompartment.NE, getRemoteFile(u5Files[0]));
		experimentU5.addReplicate(1, CellType.U, CellCompartment.C, getRemoteFile(u5Files[1]));
		proteinAccs.addAll(experimentU5.getProteinAccs());
		if (!Constants.TESTING) {
			Experiment experimentA3 = new Experiment("A2", CellType.A);
			experimentsA.add(experimentA3);
			experimentA3.addReplicate(1, CellType.A, CellCompartment.N, getRemoteFile(a2Files[0]));
			experimentA3.addReplicate(1, CellType.A, CellCompartment.NE, getRemoteFile(a2Files[1]));
			experimentA3.addReplicate(1, CellType.A, CellCompartment.C, getRemoteFile(a2Files[2]));
			proteinAccs.addAll(experimentA3.getProteinAccs());
			Experiment experimentA4 = new Experiment("A4", CellType.A);
			experimentsA.add(experimentA4);
			experimentA4.addReplicate(1, CellType.A, CellCompartment.N, getRemoteFile(a4Files[0]));
			experimentA4.addReplicate(1, CellType.A, CellCompartment.NE, getRemoteFile(a4Files[1]));
			experimentA4.addReplicate(1, CellType.A, CellCompartment.C, getRemoteFile(a4Files[2]));
			experimentA4.addReplicate(2, CellType.A, CellCompartment.N, getRemoteFile(a42Files[0]));
			experimentA4.addReplicate(2, CellType.A, CellCompartment.NE, getRemoteFile(a42Files[1]));
			experimentA4.addReplicate(2, CellType.A, CellCompartment.C, getRemoteFile(a42Files[2]));
			proteinAccs.addAll(experimentA4.getProteinAccs());

			Experiment experimentM1 = new Experiment("M1", CellType.M);
			experimentsM.add(experimentM1);
			experimentM1.addReplicate(1, CellType.M, CellCompartment.N, getRemoteFile(m1Files[0]));
			experimentM1.addReplicate(1, CellType.M, CellCompartment.NE, getRemoteFile(m1Files[1]));
			experimentM1.addReplicate(1, CellType.M, CellCompartment.C, getRemoteFile(m1Files[2]));
			experimentM1.addReplicate(2, CellType.M, CellCompartment.N, getRemoteFile(m12Files[0]));
			experimentM1.addReplicate(2, CellType.M, CellCompartment.NE, getRemoteFile(m12Files[1]));
			experimentM1.addReplicate(2, CellType.M, CellCompartment.C, getRemoteFile(m12Files[2]));
			proteinAccs.addAll(experimentM1.getProteinAccs());
			Experiment experimentM3 = new Experiment("M3", CellType.M);
			experimentsM.add(experimentM3);

			experimentM3.addReplicate(1, CellType.M, CellCompartment.N, getRemoteFile(m3Files[0]));
			experimentM3.addReplicate(1, CellType.M, CellCompartment.NE, getRemoteFile(m3Files[1]));
			experimentM3.addReplicate(1, CellType.M, CellCompartment.C, getRemoteFile(m3Files[2]));
			experimentM3.addReplicate(2, CellType.M, CellCompartment.N, getRemoteFile(m32Files[0]));
			experimentM3.addReplicate(2, CellType.M, CellCompartment.NE, getRemoteFile(m32Files[1]));
			experimentM3.addReplicate(2, CellType.M, CellCompartment.C, getRemoteFile(m32Files[2]));
			proteinAccs.addAll(experimentM3.getProteinAccs());
		}
		log.info(proteinAccs.size() + " proteins in all experiments");

		long t2 = System.currentTimeMillis();
		log.info("It took " + DatesUtil.getDescriptiveTimeFromMillisecs(t2 - t1));

	}

	private Set<String> getProteinsByCellType(CellType celltype) throws IOException {
		Set<String> ret = new HashSet<String>();
		List<Experiment> experimentList = new ArrayList<Experiment>();
		switch (celltype) {
		case A:
			experimentList.addAll(experimentsA);
			break;
		case M:
			experimentList.addAll(experimentsM);
			break;
		case U:
			experimentList.addAll(experimentsU);
			break;
		default:
			break;
		}
		for (Experiment experiment : experimentList) {
			ret.addAll(experiment.getProteinAccs());
		}

		return ret;
	}

	private RemoteSSHFileReference getRemoteFile(String remotePath) throws IOException {

		RemoteSSHFileReference ret = new RemoteSSHFileReference(hostName, userName, pass, remotefileName, null);
		ret.setRemotePath(remotePath);
		return ret;
	}

	private void run() throws IOException {
		// choose scoring function
		long t1 = System.currentTimeMillis();
		try {
			scoringFunction = new SPCScoringFunction2(this);

			// load data
			loadDatasets();

			// anotate proteins with uniprot
			// annotateProteins();

			// print scores for each cell type
			final CellType[] values = CellType.values();
			for (CellType cellType : values) {
				printScoreDistributions(cellType, false, true);
				printScoreDistributions(cellType, true, false);
			}
			// print scores for all celltypes together
			printScoreDistributions(null, false, true);
			printScoreDistributions(null, true, false);

			// compare the scores between U and A
			PairComparisonReport comparisonReportUA = compareScores(CellType.U, CellType.A);
			comparisonReportUA
					.printToFile(new File(outputFolder.getAbsolutePath() + File.separator + "U_vs_A_comparison.txt"));
			// compare the scores between U and M
			PairComparisonReport comparisonReportUM = compareScores(CellType.U, CellType.M);
			comparisonReportUM
					.printToFile(new File(outputFolder.getAbsolutePath() + File.separator + "U_vs_M_comparison.txt"));
			// compare the scores between A and M
			PairComparisonReport comparisonReportAM = compareScores(CellType.A, CellType.M);
			comparisonReportAM
					.printToFile(new File(outputFolder.getAbsolutePath() + File.separator + "A_vs_M_comparison.txt"));
			// compare the scores between U and A and M
			TripleComparisonReport comparisonReportUAM = compareScores(CellType.U, CellType.A, CellType.M);
			comparisonReportUAM.printToFile(
					new File(outputFolder.getAbsolutePath() + File.separator + "U_vs_A_vs_M_comparison.txt"));
		} finally {
			log.info("It took " + DatesUtil.getDescriptiveTimeFromMillisecs(System.currentTimeMillis() - t1));
		}
	}

	private void annotateProteins() {
		Set<String> uniprotAccs = new HashSet<String>();
		for (String rawAcc : proteinAccs) {
			String uniprotAcc = FastaParser.getUniProtACC(rawAcc);
			if (uniprotAcc == null) {
				log.warn("Not uniprot acc: " + rawAcc);
				continue;
			} else {
				uniprotAccs.add(uniprotAcc);
			}
		}
		UniprotProteinRetriever upr = Constants.upr;
		final Map<String, Protein> annotatedProteins = upr.getAnnotatedProteins(uniprotAccs);
		log.info(annotatedProteins.size() + " proteina annotated out of " + uniprotAccs.size());
	}

	private PairComparisonReport compareScores(CellType cellType1, CellType cellType2) {
		log.info(" Comparing " + cellType1 + " with " + cellType2);
		final List<Pair<String, Double>> scores1 = scoresByCellType.get(cellType1);
		Map<String, Double> scoreMap1 = getScoreMap(scores1);
		final List<Pair<String, Double>> scores2 = scoresByCellType.get(cellType2);
		Map<String, Double> scoreMap2 = getScoreMap(scores2);
		return new PairComparisonReport(this, cellType1, scoreMap1, cellType2, scoreMap2);
	}

	private TripleComparisonReport compareScores(CellType cellType1, CellType cellType2, CellType cellType3) {
		final List<Pair<String, Double>> scores1 = scoresByCellType.get(cellType1);
		Map<String, Double> scoreMap1 = getScoreMap(scores1);
		final List<Pair<String, Double>> scores2 = scoresByCellType.get(cellType2);
		Map<String, Double> scoreMap2 = getScoreMap(scores2);
		final List<Pair<String, Double>> scores3 = scoresByCellType.get(cellType3);
		Map<String, Double> scoreMap3 = getScoreMap(scores3);
		return new TripleComparisonReport(this, cellType1, scoreMap1, cellType2, scoreMap2, cellType3, scoreMap3);
	}

	private Map<String, Double> getScoreMap(List<Pair<String, Double>> scoreList) {
		Map<String, Double> ret = new HashMap<String, Double>();
		for (Pair<String, Double> pair : scoreList) {
			ret.put(pair.getFirstelement(), pair.getSecondElement());
		}
		return ret;
	}

	private void printScoreDistributions(CellType celltype, boolean peptideCounts, boolean includeRatios)
			throws IOException {
		List<Pair<String, Double>> scores = calculateScoresFromGroups(celltype);

		log.info("Printing scores for " + celltype + " to file");
		String cellTypeName = "UAM";
		if (celltype != null) {
			cellTypeName = celltype.name();
		}
		String pathname = outputFolder.getAbsolutePath() + File.separator + cellTypeName;
		if (peptideCounts) {
			pathname += "peptide_counts";
		} else {
			pathname += "spec_counts";
		}
		pathname += "_scores_distribution.txt";
		File scoreFileOutput = new File(pathname);

		FileWriter fw = null;
		try {
			fw = new FileWriter(scoreFileOutput);
			writeHeaders(fw, celltype, peptideCounts, includeRatios);
			int num = 1;
			int previousPercentaje = 0;
			for (Pair<String, Double> pair : scores) {
				int percentage = Double.valueOf(num * 1.0 / scores.size() * 100.0).intValue();
				if (percentage != previousPercentaje) {
					previousPercentaje = percentage;
					log.info(previousPercentaje + "% proteins printed");
				}

				if (!Double.isNaN(pair.getSecondElement())) {
					String rawAcc = pair.getFirstelement();
					if (rawAcc.contains("P29788,Q8K3F2[CONCLUSIVE]")) {
						log.info(rawAcc);
					}
					if (rawAcc.contains("[")) {
						rawAcc = rawAcc.substring(0, rawAcc.indexOf("["));
					}
					String rawAccString = getAccessionStringByEvidence(rawAcc);
					String geneNameString = getGeneNameString(rawAcc);
					if (geneNameString.toLowerCase().contains("nab2")) {
						System.out.println(pair);
					}
					String proteinNameString = getProteinNameString(rawAcc);
					// testing
					if (ControlNE.isControl(rawAcc)) {
						double score = scoringFunction.getScore(getGroup(rawAcc), celltype);
						if (score < 1) {
							score = scoringFunction.getScore(getGroup(rawAcc), celltype);
						}
					}
					// end testing
					String valid = isValid(rawAcc) ? "VALID" : "FILTERED";
					fw.write(num++ + "\t" + valid + "\t" + rawAccString + "\t" + geneNameString + "\t"
							+ proteinNameString + "\t" + getProteinEvidences(rawAcc) + "\t" + pair.getSecondElement()
							+ "\t" + ControlNE.isControl(rawAcc) + "\t" + Constants.isDecoy(getAccs(rawAcc)) + "\t"
							+ getTransmembraneRegion(rawAcc) + "\t");
					for (Experiment experiment : getAllExperiments()) {
						if (celltype == null || experiment.getCellType() == celltype) {
							// print the averages
							for (CellCompartment cellCompartment : CellCompartment.values()) {

								if (peptideCounts) {
									fw.write(experiment.getAvgPeptideCount(getGroup(rawAcc), cellCompartment, true)
											+ "\t");
								} else {
									fw.write(experiment.getAvgSpectralCount(getGroup(rawAcc), cellCompartment, true)
											+ "\t");
								}

								if (cellCompartment != Constants.cellCompartmentToStudy) {
									if (includeRatios) {
										double spcRatio = experiment.getSPCRatio(getGroup(rawAcc),
												Constants.cellCompartmentToStudy, cellCompartment, true);

										fw.write(Constants.formatInfinity(ScoringFunction.getLog2Ratio(spcRatio))
												+ "\t");
									}
								}
							}

							fw.write(experiment.printColumnsForProteinGroup(getGroup(rawAcc), peptideCounts));

						}
					}

					fw.write("\n");

				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			log.error(e.getMessage());
		} finally {
			if (fw != null) {
				fw.close();
			}
			log.info("File created at " + scoreFileOutput.getAbsolutePath());
		}

	}

	private boolean getTransmembraneRegion(String rawAcc) throws IOException {

		List<Boolean> validArray = filterAccessionsByEvidence(rawAcc);
		int index = 0;
		for (String acc : getAccs(rawAcc)) {
			if (!validArray.get(index++)) {
				continue;
			}
			final Map<String, Protein> annotatedProtein = Constants.upr.getAnnotatedProtein(acc);
			if (annotatedProtein.containsKey(acc)) {
				final Protein protein2 = annotatedProtein.get(acc);
				if (protein2 != null) {
					Set<ProteinAnnotation> annotations = protein2.getAnnotations();
					for (ProteinAnnotation proteinAnnotation : annotations) {
						if (proteinAnnotation.getAnnotationType().getKey().equals("transmembrane region")) {
							return true;
						}
					}

				}
			}
		}
		return false;
	}

	private String getAccessionStringByEvidence(String rawAcc) throws IOException {
		List<Boolean> validArray = filterAccessionsByEvidence(rawAcc);
		int index = 0;
		StringBuilder sb = new StringBuilder();
		for (String acc : getAccs(rawAcc)) {
			if (!validArray.get(index++)) {
				continue;
			}
			if (!"".equals(sb.toString())) {
				sb.append(Constants.SEPARATOR);
			}
			sb.append(acc);

		}
		return sb.toString();
	}

	private List<Boolean> filterAccessionsByEvidence(String rawAcc) throws IOException {

		Set<String> accs = getAccs(rawAcc);
		final Map<String, Entry> uniprotEntries = Constants.upr.getAnnotatedUniprotEntries(accs);
		ProteinGroup group = getGroup(rawAcc);
		List<Boolean> ret = new ArrayList<Boolean>();
		List<Boolean> groupEvidenceArray = new ArrayList<Boolean>();
		List<Boolean> uniprotEvidenceArray = new ArrayList<Boolean>();
		// only swissprot is valid
		boolean thereIsASwissProt = false;
		boolean thereIsAConclusiveProt = false;
		for (String acc : accs) {
			boolean valid = false;
			try {

				ProteinEvidence evidence = getEvidence(group, acc);
				if (evidence == ProteinEvidence.CONCLUSIVE) {
					groupEvidenceArray.add(true);
					thereIsAConclusiveProt = true;

					valid = true;
				} else {
					groupEvidenceArray.add(false);
					if (uniprotEntries.containsKey(acc)) {
						final Entry protein2 = uniprotEntries.get(acc);
						if (protein2 != null) {
							String dataset = protein2.getDataset();
							if (dataset != null) {
								if (dataset.toLowerCase().equals("swiss-prot")) {
									thereIsASwissProt = true;
									uniprotEvidenceArray.add(true);
									valid = true;
								} else {
									uniprotEvidenceArray.add(false);
								}
							}
						} else {
							uniprotEvidenceArray.add(false);
						}
					} else {
						uniprotEvidenceArray.add(false);
					}
				}
			} finally {
				ret.add(valid);
			}
		}
		boolean allSwissprot = true;
		for (Boolean uniprotEvidence : uniprotEvidenceArray) {
			if (!uniprotEvidence) {
				allSwissprot = false;
			}
		}

		if (!thereIsASwissProt && !thereIsAConclusiveProt) {
			ret.set(0, true);
		}
		if (allSwissprot) {
			if (groupEvidenceArray.size() != accs.size()) {
				log.info(rawAcc);
			}
			// only report the conclusive ones
			return groupEvidenceArray;
		}
		if (ret.size() != accs.size()) {
			log.info(rawAcc);
		}
		return ret;
	}

	private ProteinEvidence getEvidence(ProteinGroup group, String acc) {
		for (GroupableProtein groupableProtein : group) {
			if (groupableProtein.getAccession().equals(acc)) {
				return groupableProtein.getEvidence();
			}
		}
		return null;
	}

	private void writeHeaders(FileWriter fw, CellType celltype, boolean peptideCounts, boolean includeRatios)
			throws IOException {
		String dataType = "SPC";
		if (peptideCounts) {
			dataType = "PEPC";
		}
		// write the header
		fw.write(
				"NUM\tVALID\tACC\tGene\tprotein description\tPROTEIN_EVIDENCE\tSCORE\tKnown NE\tdecoy\tTransmembrane region\t");
		List<Experiment> experimentList = new ArrayList<Experiment>();
		experimentList.addAll(experimentsU);
		experimentList.addAll(experimentsA);
		experimentList.addAll(experimentsM);
		for (Experiment experiment : experimentList) {
			if (celltype == null || experiment.getCellType() == celltype) {
				// print the averages
				for (CellCompartment cellCompartment : CellCompartment.values()) {
					fw.write(experiment.getName() + cellCompartment.name() + "_AVG_" + dataType + "\t");
					if (includeRatios) {
						if (cellCompartment != Constants.cellCompartmentToStudy) {
							fw.write(experiment.getName() + Constants.cellCompartmentToStudy + "/"
									+ experiment.getName() + cellCompartment + "\t");
						}
					}
				}
				fw.write(experiment.printHeader());
			}
		}

		fw.write("\n");

	}

	public Set<String> getAccs(String rawAcc) {
		Set<String> accs = new HashSet<String>();

		// remove the [evidence] at the end
		if (rawAcc.contains("[")) {
			rawAcc = rawAcc.substring(0, rawAcc.indexOf("["));
		}
		if (rawAcc.contains(",")) {
			final String[] split = rawAcc.split(",");
			for (String acc : split) {
				if (acc.contains("|")) {
					log.info(rawAcc);
				}
				accs.add(acc);
			}
		} else {
			accs.add(rawAcc);
		}
		return accs;

	}

	public String getProteinNameString(String rawAcc) throws IOException {

		String proteinName = "";
		List<Boolean> validArray = filterAccessionsByEvidence(rawAcc);
		int index = 0;
		for (String acc : getAccs(rawAcc)) {
			if (!validArray.get(index++)) {
				continue;
			}
			final Map<String, Protein> annotatedProtein = Constants.upr.getAnnotatedProtein(acc);
			if (annotatedProtein.containsKey(acc)) {
				final Protein protein2 = annotatedProtein.get(acc);
				if (protein2 != null) {
					if (protein2.getPrimaryAccession() != null
							&& protein2.getPrimaryAccession().getDescription() != null) {
						if (!"".equals(proteinName)) {
							proteinName += Constants.SEPARATOR;
						}
						proteinName += protein2.getPrimaryAccession().getDescription();
					}

				}
			}
		}
		return proteinName;
	}

	public String getGeneNameString(String rawAcc) throws IOException {

		List<String> geneNames = new ArrayList<String>();
		List<Boolean> validArray = filterAccessionsByEvidence(rawAcc);
		int index = 0;
		for (String acc : getAccs(rawAcc)) {
			if (!validArray.get(index++)) {
				continue;
			}
			final Map<String, Protein> annotatedProtein = Constants.upr.getAnnotatedProtein(acc);
			if (annotatedProtein.containsKey(acc)) {
				final Protein protein2 = annotatedProtein.get(acc);
				if (protein2 != null) {
					boolean onePrimary = false;
					if (protein2.getGenes() != null) {
						for (Gene gene : protein2.getGenes()) {
							if (gene.getGeneType().equals("primary")) {
								onePrimary = true;
								if (!geneNames.contains(gene.getGeneID())) {
									geneNames.add(gene.getGeneID());
								}
							}
						}
						if (!onePrimary) {
							if (protein2.getGenes() != null && !protein2.getGenes().isEmpty()) {
								geneNames.add(protein2.getGenes().iterator().next().getGeneID());
							}
						}
					}
				}
			}
		}
		return getGeneName(geneNames);
	}

	private List<Experiment> getAllExperiments() {
		List<Experiment> list = new ArrayList<Experiment>();
		list.addAll(experimentsA);
		list.addAll(experimentsM);
		list.addAll(experimentsU);
		return list;
	}

	private ProteinGroup getGroup(String groupKey) throws IOException {
		List<ProteinGroup> proteinGroups = getProteinGroups().stream()
				.filter(group -> group.getKey().contains(groupKey) && !group.getKey().toLowerCase().contains("reverse"))
				.collect(Collectors.toList());
		if (proteinGroups.size() > 1) {
			log.warn("More than one group with key " + groupKey);
		}
		if (proteinGroups.isEmpty()) {
			proteinGroups = getProteinGroups().stream().filter(group -> group.getKey().contains(groupKey))
					.collect(Collectors.toList());
			if (!proteinGroups.isEmpty()) {
				return proteinGroups.get(0);
			}
		}
		if (!proteinGroups.isEmpty()) {
			return proteinGroups.get(0);
		}
		return null;
	}

	private String getProteinEvidences(String groupKey) throws IOException {
		List<String> proteinAccs = new ArrayList<String>();
		if (groupKey == null) {
			return "";
		}
		if (groupKey.contains("[")) {
			groupKey = groupKey.substring(0, groupKey.indexOf("["));
		}
		if (groupKey.contains(",")) {
			String[] split = groupKey.split(",");
			for (String string : split) {
				proteinAccs.add(string);
			}
		} else {
			proteinAccs.add(groupKey);
		}

		ProteinGroup group = getGroup(groupKey);
		StringBuilder sb = new StringBuilder();
		if (group != null) {
			for (String proteinAcc : proteinAccs) {
				ProteinEvidence evidence = null;
				for (GroupableProtein groupableProtein : group) {
					if (groupableProtein.getAccession().equals(proteinAcc)) {
						evidence = groupableProtein.getEvidence();
					}
				}
				if (!"".equals(sb.toString())) {
					sb.append(Constants.SEPARATOR);
				}
				if (evidence != null) {
					sb.append(evidence);
				} else {
					log.info(" I didnt find " + proteinAcc + "  in group  " + group);
					for (GroupableProtein groupableProtein : group) {
						log.info(groupableProtein.getAccession());
					}
					sb.append("");
				}
			}

		}
		return sb.toString();
	}

	private String getGeneName(List<String> geneNames) {
		StringBuilder sb = new StringBuilder();
		Set<String> set = new HashSet<String>();
		for (String geneName : geneNames) {
			if (set.contains(geneName)) {
				continue;
			}
			if (!"".equals(sb.toString())) {
				sb.append(Constants.SEPARATOR);
			}
			sb.append(geneName);
			set.add(geneName);
		}
		return sb.toString();
	}

	/**
	 * @return the experimentU
	 */
	public List<Experiment> getExperimentsU() {
		return experimentsU;
	}

	/**
	 * @return the experimentA
	 */
	public List<Experiment> getExperimentsA() {
		return experimentsA;
	}

	/**
	 * @return the experimentM
	 */
	public List<Experiment> getExperimentsM() {
		return experimentsM;
	}

	private Filter getGOFilter() {

		if (Constants.GO_FILTER && goFilter == null) {
			goFilter = new GOFilter("ASDF");

		}

		return goFilter;
	}

	public boolean isValid(String rawAcc) throws IOException {
		final Filter goFilter2 = getGOFilter();

		final ProteinGroup group = getGroup(rawAcc);
		if (group != null) {
			for (GroupableProtein groupableProtein : group) {
				if (goFilter2 != null && !goFilter2.isValid((Protein) groupableProtein)) {
					return false;
				}

			}
			return true;
		}

		return true;
	}
}
