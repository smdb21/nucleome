package edu.scripps.yates.nucleome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.time.DateFormatUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinRetriever;
import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.nucleome.filters.Filter;
import edu.scripps.yates.nucleome.filters.GOFilter;
import edu.scripps.yates.nucleome.filters.KeratinFilter;
import edu.scripps.yates.nucleome.model.CellCompartment;
import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.nucleome.model.ControlNE;
import edu.scripps.yates.nucleome.model.Experiment;
import edu.scripps.yates.nucleome.model.Fractionation;
import edu.scripps.yates.nucleome.model.Replicate;
import edu.scripps.yates.utilities.dates.DatesUtil;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.grouping.GroupablePSM;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.PAnalyzer;
import edu.scripps.yates.utilities.grouping.ProteinEvidence;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.proteomicsmodel.Gene;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.ProteinAnnotation;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.set.hash.THashSet;

public class _4DNucleomeAnalyzer {
	private final static Logger log = Logger.getLogger(_4DNucleomeAnalyzer.class);

	public static void main(String[] args) {
		_4DNucleomeAnalyzer analyzer;
		try {

			analyzer = new _4DNucleomeAnalyzer();

			////////////////////////////////////////////////////////////
			// PARAMETERS
			Constants.includeNegativeScoring = false;
			Constants.MIN_PEPTIDES_PER_PROTEIN = 1;
			Constants.MIN_PSM_PER_PROTEIN = 2;
			Constants.MIN_AVG_SPC = 3;
			Constants.GO_FILTER = false;
			Constants.cellCompartmentToStudy = CellCompartment.NE;
			Constants.TESTING = false;
			Constants.ENRICHMENT_SCORE_THRESHOLD = null;// 3.0;
			Constants.DATASET_PATHS_FILE = datasetsPathsFile;
			Constants.MIN_TOTAL_SPC = 5;
			Constants.MAX_TEST_PROTEINS = 200000;
			Constants.writeCombinedDistribution = false;// UAM
			Constants.compareScores = false;
			UniprotProteinRetriever.enableCache = true;
			// scoringFunction = new ScoringFunctionByNE_SPC_Percentage(this);
			scoringFunction = new ScoringFunctionByNE_NSAF_Percentage(analyzer);
			////////////////////////////////////////////////////////////

			analyzer.run();

			System.out.println("DONE");
			System.exit(0);
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("ERROR: " + e.getMessage());
		}
		System.exit(-1);
	}

	private final List<Experiment> experimentsU = new ArrayList<Experiment>();
	private final List<Experiment> experimentsA = new ArrayList<Experiment>();
	private final List<Experiment> experimentsM = new ArrayList<Experiment>();
	private static ScoringFunction scoringFunction;
	private final String hostName = "jaina.scripps.edu";;
	private final String userName = "salvador";
	private final String pass = "Natjeija21";
	private final String remotefileName = "DTASelect-filter.txt";
	private static final String datasetsPathsFile = "z:\\share\\Salva\\data\\4D_Nucleome\\SwissProt_1FDR\\SwissProt_1FDR_data_paths.txt";
	private static final String datasetsPhosphoPathsFile = "z:\\share\\Salva\\data\\4D_Nucleome\\datasets_paths_phospho.txt";
	private final File outputFolder = new File(new File(datasetsPathsFile).getParent() + File.separator + "output");

	private final Map<CellType, List<Pair<String, Double>>> scoresByCellType = new HashMap<CellType, List<Pair<String, Double>>>();
	private GOFilter goFilter;
	private KeratinFilter keratinFilter;
	private int fileNum = 1;
	private List<ProteinGroup> proteinGroups;
	private HashSet<String> proteinAccs;
	private Map<String, Integer> totalSPCs = new HashMap<String, Integer>();
	private Map<String, ProteinGroup> groupsByRawAcc = new HashMap<String, ProteinGroup>();
	private Set<GroupableProtein> groupableProteins = new THashSet<GroupableProtein>();

	public _4DNucleomeAnalyzer() throws IOException {

	}

	/**
	 * Returns a list of pairs (Protein ACC - Score), sorted by the score
	 *
	 * @return
	 * @throws IOException
	 */
	private List<Pair<String, Double>> calculateScoresFromGroups(CellType celltype) throws IOException {
		log.info("Calculating scores for " + getAllAccs(celltype).size() + " proteins in " + celltype);
		List<Pair<String, Double>> scores = new ArrayList<Pair<String, Double>>();
		List<ProteinGroup> proteinGroups = getProteinGroups(celltype);
		for (ProteinGroup proteinGroup : proteinGroups) {
			if (proteinGroup.getEvidence() == ProteinEvidence.NONCONCLUSIVE) {
				continue;
			}
			String rawAcc = proteinGroup.getKey();
			List<String> filteredAcessions = new ArrayList<String>();
			if (rawAcc.contains("P05213")) {
				log.info(rawAcc);
				getAccessionStringByEvidence(rawAcc, proteinGroup, celltype);
			}
			String filteredAcc = getAccessionStringByEvidence(rawAcc, proteinGroup, celltype);
			if (filteredAcc.contains(",")) {
				String[] split = filteredAcc.split(",");
				for (String acc : split) {
					filteredAcessions.add(acc);
				}
			} else {
				filteredAcessions.add(filteredAcc);
			}
			double score = scoringFunction.getScore(filteredAcessions, celltype);

			Pair<String, Double> pair = new Pair<String, Double>(filteredAcc, score);
			scores.add(pair);

		}
		log.info("Sorting scores...");
		// sort them by the score
		Collections.sort(scores, new Comparator<Pair<String, Double>>() {

			@Override
			public int compare(Pair<String, Double> o1, Pair<String, Double> o2) {

				try {
					String rawAcc1 = o1.getFirstelement();
					if (rawAcc1.contains("[")) {
						rawAcc1 = rawAcc1.substring(0, rawAcc1.indexOf("["));
					}
					if (rawAcc1.contains("P05213")) {
						log.info(rawAcc1);
					}
					double score1 = o1.getSecondElement();
					int totalSPC1 = getTotalSPC(rawAcc1, celltype);
					boolean valid1 = isValid(rawAcc1, totalSPC1) ? true : false;
					//
					String rawAcc2 = o2.getFirstelement();
					if (rawAcc2.contains("[")) {
						rawAcc2 = rawAcc2.substring(0, rawAcc2.indexOf("["));
					}
					double score2 = o2.getSecondElement();
					int totalSPC2 = getTotalSPC(rawAcc2, celltype);
					boolean valid2 = isValid(rawAcc2, totalSPC2) ? true : false;
					//

					if (valid1 && !valid2) {
						return -1;
					} else if (!valid1 && valid2) {
						return 1;
					} else {
						int scoreComparison = Double.compare(score2, score1);
						if (scoreComparison == 0) {
							// by total SPC
							return Integer.compare(totalSPC2, totalSPC1);
						} else {
							return scoreComparison;
						}
					}

				} catch (IOException e) {
					return 0;
				}
			}
		});
		log.info(scores.size() + " scores calculated");

		scoresByCellType.put(celltype, scores);
		return scores;
	}

	private List<ProteinGroup> getProteinGroups(CellType cellType) throws IOException {
		if (proteinGroups == null) {
			PAnalyzer panalyzer = new PAnalyzer(true);
			Set<GroupableProtein> groupableProteins = getAllGroupableProteins(cellType);

			log.info("Grouping " + groupableProteins.size() + " proteins from " + cellType);
			List<ProteinGroup> proteinGroupsTMP = panalyzer.run(groupableProteins);
			log.info(proteinGroupsTMP.size() + " protein groups");
			proteinGroups = new ArrayList<ProteinGroup>();
			for (ProteinGroup proteinGroup : proteinGroupsTMP) {
				if (proteinGroup.getEvidence() != ProteinEvidence.NONCONCLUSIVE) {
					proteinGroups.add(proteinGroup);
				}
			}
			log.info(proteinGroups.size() + " protein groups after removing NonConclusive proteins");
		}
		return proteinGroups;
	}

	private Set<GroupableProtein> getAllGroupableProteins(CellType cellType) throws IOException {
		if (groupableProteins.isEmpty()) {
			for (Experiment experiment : getExperimentsA()) {
				if (cellType != null && !experiment.getCellType().equals(cellType)) {
					continue;
				}
				for (Protein protein : experiment.getProteins()) {
					groupableProteins.add(protein);
				}
			}
			for (Experiment experiment : getExperimentsM()) {
				if (cellType != null && !experiment.getCellType().equals(cellType)) {
					continue;
				}
				for (Protein protein : experiment.getProteins()) {
					groupableProteins.add(protein);
				}
			}
			for (Experiment experiment : getExperimentsU()) {
				if (cellType != null && !experiment.getCellType().equals(cellType)) {
					continue;
				}
				for (Protein protein : experiment.getProteins()) {
					groupableProteins.add(protein);
				}
			}
		}
		return groupableProteins;
	}

	private void loadDatasets() throws IOException {
		log.info("Loading datasets");
		long t1 = System.currentTimeMillis();

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
		// String[] u5Files = dataPaths.getFiles("U5");
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

		Experiment experimentU3 = new Experiment("U3", CellType.U);
		experimentsU.add(experimentU3);
		experimentU3.addReplicate(1, CellType.U, CellCompartment.N, getRemoteFile(u3Files[0]));
		experimentU3.addReplicate(1, CellType.U, CellCompartment.NE, getRemoteFile(u3Files[1]));
		experimentU3.addReplicate(1, CellType.U, CellCompartment.C, getRemoteFile(u3Files[2]));

		Experiment experimentU4 = new Experiment("U4", CellType.U);
		experimentsU.add(experimentU4);
		experimentU4.addReplicate(1, CellType.U, CellCompartment.N, getRemoteFile(u4Files[0]));
		experimentU4.addReplicate(1, CellType.U, CellCompartment.NE, getRemoteFile(u4Files[1]));
		experimentU4.addReplicate(1, CellType.U, CellCompartment.C, getRemoteFile(u4Files[2]));
		experimentU4.addReplicate(2, CellType.U, CellCompartment.N, getRemoteFile(u42Files[0]));
		experimentU4.addReplicate(2, CellType.U, CellCompartment.NE, getRemoteFile(u42Files[1]));
		experimentU4.addReplicate(2, CellType.U, CellCompartment.C, getRemoteFile(u42Files[2]));

		// Experiment experimentU5 = new Experiment("U5", CellType.U);
		// experimentsU.add(experimentU5);
		// experimentU5.addReplicate(1, CellType.U, CellCompartment.NE,
		// getRemoteFile(u5Files[0]));
		// experimentU5.addReplicate(1, CellType.U, CellCompartment.C,
		// getRemoteFile(u5Files[1]));

		if (!Constants.TESTING) {
			Experiment experimentA3 = new Experiment("A2", CellType.A);
			experimentsA.add(experimentA3);
			experimentA3.addReplicate(1, CellType.A, CellCompartment.N, getRemoteFile(a2Files[0]));
			experimentA3.addReplicate(1, CellType.A, CellCompartment.NE, getRemoteFile(a2Files[1]));
			experimentA3.addReplicate(1, CellType.A, CellCompartment.C, getRemoteFile(a2Files[2]));

			Experiment experimentA4 = new Experiment("A4", CellType.A);
			experimentsA.add(experimentA4);
			experimentA4.addReplicate(1, CellType.A, CellCompartment.N, getRemoteFile(a4Files[0]));
			experimentA4.addReplicate(1, CellType.A, CellCompartment.NE, getRemoteFile(a4Files[1]));
			experimentA4.addReplicate(1, CellType.A, CellCompartment.C, getRemoteFile(a4Files[2]));
			experimentA4.addReplicate(2, CellType.A, CellCompartment.N, getRemoteFile(a42Files[0]));
			experimentA4.addReplicate(2, CellType.A, CellCompartment.NE, getRemoteFile(a42Files[1]));
			experimentA4.addReplicate(2, CellType.A, CellCompartment.C, getRemoteFile(a42Files[2]));

			Experiment experimentM1 = new Experiment("M1", CellType.M);
			experimentsM.add(experimentM1);
			experimentM1.addReplicate(1, CellType.M, CellCompartment.N, getRemoteFile(m1Files[0]));
			experimentM1.addReplicate(1, CellType.M, CellCompartment.NE, getRemoteFile(m1Files[1]));
			experimentM1.addReplicate(1, CellType.M, CellCompartment.C, getRemoteFile(m1Files[2]));
			experimentM1.addReplicate(2, CellType.M, CellCompartment.N, getRemoteFile(m12Files[0]));
			experimentM1.addReplicate(2, CellType.M, CellCompartment.NE, getRemoteFile(m12Files[1]));
			experimentM1.addReplicate(2, CellType.M, CellCompartment.C, getRemoteFile(m12Files[2]));

			Experiment experimentM3 = new Experiment("M3", CellType.M);
			experimentsM.add(experimentM3);

			experimentM3.addReplicate(1, CellType.M, CellCompartment.N, getRemoteFile(m3Files[0]));
			experimentM3.addReplicate(1, CellType.M, CellCompartment.NE, getRemoteFile(m3Files[1]));
			experimentM3.addReplicate(1, CellType.M, CellCompartment.C, getRemoteFile(m3Files[2]));
			experimentM3.addReplicate(2, CellType.M, CellCompartment.N, getRemoteFile(m32Files[0]));
			experimentM3.addReplicate(2, CellType.M, CellCompartment.NE, getRemoteFile(m32Files[1]));
			experimentM3.addReplicate(2, CellType.M, CellCompartment.C, getRemoteFile(m32Files[2]));

		}

		long t2 = System.currentTimeMillis();
		log.info("It took " + DatesUtil.getDescriptiveTimeFromMillisecs(t2 - t1));

	}

	private File getRemoteFile(String remotePath) throws IOException {
		File localFile = new File(outputFolder + File.separator + "data_files" + File.separator
				+ FilenameUtils.getBaseName(remotePath + remotefileName) + "_" + fileNum++ + "."
				+ FilenameUtils.getExtension(remotefileName));
		if (localFile.exists()) {
			return localFile;
		}
		RemoteSSHFileReference ret = new RemoteSSHFileReference(hostName, userName, pass, remotefileName, null);
		ret.setRemotePath(remotePath);

		ret.setOutputFile(localFile);
		return ret.getRemoteFile();
	}

	public void run() throws IOException {
		// choose scoring function
		long t1 = System.currentTimeMillis();
		try {
			// load data
			loadDatasets();
			// print scores for each cell type
			final CellType[] values = CellType.values();
			for (CellType cellType : values) {
				// restart field variables
				proteinAccs = null;
				proteinGroups = null;
				totalSPCs.clear();
				groupsByRawAcc.clear();
				groupableProteins.clear();

				// annotate proteins with uniprot
				annotateProteins(cellType);
				writeScoreDistributions(cellType);

				// writeScoreDistributions(cellType, DataType.NSAF, true &&
				// Constants.printRatios);
				// writeScoreDistributions(cellType, DataType.PEPC, false &&
				// Constants.printRatios);

			}
			// print scores for all celltypes together
			// writeScoreDistributions(null, DataType.NSAF, true &&
			// Constants.printRatios);
			if (Constants.writeCombinedDistribution) {
				writeScoreDistributions(null);
			}
			// writeScoreDistributions(null, DataType.PEPC, false &&
			// Constants.printRatios);

			if (Constants.compareScores) {
				// compare the scores between U and A
				PairComparisonReport comparisonReportUA = compareScores(CellType.U, CellType.A);
				comparisonReportUA.printToFile(
						new File(outputFolder.getAbsolutePath() + File.separator + "U_vs_A_comparison.txt"));
				// compare the scores between U and M
				PairComparisonReport comparisonReportUM = compareScores(CellType.U, CellType.M);
				comparisonReportUM.printToFile(
						new File(outputFolder.getAbsolutePath() + File.separator + "U_vs_M_comparison.txt"));
				// compare the scores between A and M
				PairComparisonReport comparisonReportAM = compareScores(CellType.A, CellType.M);
				comparisonReportAM.printToFile(
						new File(outputFolder.getAbsolutePath() + File.separator + "A_vs_M_comparison.txt"));
				// compare the scores between U and A and M
				TripleComparisonReport comparisonReportUAM = compareScores(CellType.U, CellType.A, CellType.M);
				comparisonReportUAM.printToFile(
						new File(outputFolder.getAbsolutePath() + File.separator + "U_vs_A_vs_M_comparison.txt"));
			}
		} finally {
			log.info("It took " + DatesUtil.getDescriptiveTimeFromMillisecs(System.currentTimeMillis() - t1));
		}
	}

	private void annotateProteins(CellType cellType) throws IOException {
		Set<String> uniprotAccs = new HashSet<String>();
		for (String rawAcc : getAllAccs(cellType)) {
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
		log.info(annotatedProteins.size() + " proteins annotated out of " + uniprotAccs.size());
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

	private void writeScoreDistributions(CellType celltype) throws IOException {
		if (!Constants.printScoreDistributions) {
			return;
		}
		List<Pair<String, Double>> scores = calculateScoresFromGroups(celltype);

		log.info("Printing scores for " + celltype + " to file");
		String cellTypeName = "UAM";
		if (celltype != null) {
			cellTypeName = celltype.name();
		}
		String formatedDate = DateFormatUtils.format(new Date(), "yyyy-MM-dd_HH-mm");
		String pathname = outputFolder.getAbsolutePath() + File.separator + formatedDate + "_" + cellTypeName
				+ "_scores_distribution.txt";

		File scoreFileOutput = new File(pathname);
		if (!scoreFileOutput.getParentFile().exists()) {
			log.info("Creating file " + scoreFileOutput.getAbsolutePath());
			scoreFileOutput.getParentFile().mkdirs();
		}
		if (scoreFileOutput.exists()) {
			log.info("Overriding file " + scoreFileOutput.getAbsolutePath());
		}
		FileWriter fw = null;
		try {
			ProgressCounter counter = new ProgressCounter(scores.size(), ProgressPrintingType.PERCENTAGE_STEPS, 0);
			fw = new FileWriter(scoreFileOutput);
			writeHeaders(fw, celltype);
			int num = 1;
			for (Pair<String, Double> pair : scores) {
				counter.increment();
				String percentage = counter.printIfNecessary();
				if (!"".equals(percentage)) {
					log.info(percentage);
				}
				if (Constants.TESTING && counter.getCount() == Constants.MAX_TEST_PROTEINS) {
					break;
				}
				if (!Double.isNaN(pair.getSecondElement())) {
					String rawAcc = pair.getFirstelement();
					if (rawAcc.contains("[")) {
						rawAcc = rawAcc.substring(0, rawAcc.indexOf("["));
					}
					String rawAccString = getAccessionStringByEvidence(rawAcc, null, celltype);
					if (rawAccString.equals("Q91W53")) {
						log.info(rawAcc);
					}
					List<String> filteredAccessions = getAccs(rawAccString);
					String geneNameString = null;
					// TODO
					// if (false) {
					geneNameString = getGeneNameString(rawAcc, celltype);
					// }

					String proteinNameString = getProteinNameString(rawAcc, celltype);
					double score = pair.getSecondElement();
					// end testing
					int totalSPC = getTotalSPC(rawAcc, celltype);
					String valid = isValid(rawAcc, totalSPC) ? "VALID" : "FILTERED";

					fw.write(num++ + "\t" + valid + "\t" + rawAccString + "\t" + geneNameString + "\t"
							+ proteinNameString + "\t" + getFilteredProteinEvidences(rawAcc, celltype) + "\t" + score
							+ "\t" + ControlNE.isControl(rawAcc) + "\t" + getTransmembraneRegion(rawAcc, celltype)
							+ "\t" + totalSPC + "\t");
					if (celltype != null) {
						for (Experiment experiment : getAllExperiments()) {
							if (experiment.getCellType() == celltype) {
								// print the averages
								for (CellCompartment cellCompartment : CellCompartment.values()) {
									List<Replicate> replicates = experiment.getSortedReplicates();
									// SPCs
									for (Replicate replicate : replicates) {
										Fractionation fractionation = replicate.getFractionation(cellCompartment);
										if (fractionation != null) {
											int repSPC = fractionation.getSpectralCount(filteredAccessions, true);
											fw.write(repSPC + "\t");
										}
									}
								}
							}
						}
						for (Experiment experiment : getAllExperiments()) {
							if (experiment.getCellType() == celltype) {
								// print the averages
								for (CellCompartment cellCompartment : CellCompartment.values()) {
									List<Replicate> replicates = experiment.getSortedReplicates();
									// NSAFs
									for (Replicate replicate : replicates) {
										Fractionation fractionation = replicate.getFractionation(cellCompartment);
										if (fractionation != null) {
											double repNSAF = fractionation.getAverageNSAF(filteredAccessions, true);
											fw.write(repNSAF + "\t");
										}
									}
								}
							}
						}
					} else {
						throw new IllegalArgumentException("cell type null is not supported");
					}
					fw.write("\n");

				}
			}
		} catch (

		Exception e)

		{
			e.printStackTrace();
			log.error(e.getMessage());
		} finally

		{
			if (fw != null) {
				fw.close();
			}
			log.info("File created at " + scoreFileOutput.getAbsolutePath());
		}

	}

	public int getTotalSPC(String rawAcc, CellType cellType) throws IOException {
		if (!totalSPCs.containsKey(rawAcc)) {
			Set<String> individualAccs = new HashSet<String>();
			if (rawAcc.contains(",")) {
				// because if this has a comma it is because
				// it is an indistinguisable group
				String[] split = rawAcc.split(",");
				for (String string : split) {
					individualAccs.add(string);
				}
			} else {
				individualAccs.add(rawAcc);
			}
			Set<GroupablePSM> psms = new HashSet<GroupablePSM>();

			Set<GroupableProtein> proteins = getAllGroupableProteins(cellType);
			for (GroupableProtein groupableProtein : proteins) {
				if (rawAcc.contains(groupableProtein.getAccession())) {
					psms.addAll(groupableProtein.getGroupablePSMs());
				}
			}
			int numPSMs = psms.size();
			int numPSMsTMP = 0;
			for (Experiment experiment : getAllExperiments()) {
				if (experiment.getCellType() == cellType) {
					// print the averages
					for (CellCompartment cellCompartment : CellCompartment.values()) {
						List<Replicate> replicates = experiment.getSortedReplicates();
						// SPCs
						for (Replicate replicate : replicates) {
							Fractionation fractionation = replicate.getFractionation(cellCompartment);
							if (fractionation != null) {
								for (String acc : individualAccs) {

									int spectralCount = fractionation.getSpectralCount(acc, true);
									numPSMsTMP += spectralCount;
									// log.info(acc + " " + spectralCount + " in
									// " + fractionation.getName());
									if (spectralCount > 0) {
										// we check all of them but we only
										// count one
										// in some dtaselects may appear the
										// first one and in other may appear the
										// second,
										// so we have to check all of them
										break;
									}
								}
							}
						}
					}
				}
			}
			if (numPSMs != numPSMsTMP) {
				log.warn("cuidado ");
			}
			totalSPCs.put(rawAcc, numPSMs);
		}
		return totalSPCs.get(rawAcc);
	}

	private boolean getTransmembraneRegion(String rawAcc, CellType cellType) throws IOException {

		List<Boolean> validArray = filterAccessionsByEvidence(rawAcc, null, cellType);
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

	private String getAccessionStringByEvidence(String rawAcc, ProteinGroup proteinGroup, CellType cellType)
			throws IOException {
		List<Boolean> validArray = filterAccessionsByEvidence(rawAcc, proteinGroup, cellType);
		int index = 0;
		StringBuilder sb = new StringBuilder();
		for (String acc : getAccs(rawAcc)) {
			if (!validArray.get(index++)) {
				continue;
			}
			if (!"".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(acc);

		}
		return sb.toString();
	}

	private List<Boolean> filterAccessionsByEvidence(String rawAcc, ProteinGroup group, CellType cellType)
			throws IOException {

		List<String> accs = getAccs(rawAcc);
		final Map<String, Entry> uniprotEntries = Constants.upr.getAnnotatedUniprotEntries(accs);
		if (group == null) {
			group = getGroup(rawAcc, cellType);
		}
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
				} else if (evidence == ProteinEvidence.NONCONCLUSIVE) {
					groupEvidenceArray.add(false);
					valid = false;
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
			if (thereIsAConclusiveProt) {
				// only report the conclusive ones
				return groupEvidenceArray;
			} else {
				// do not report nonconclusive

				if (uniprotEvidenceArray.size() == ret.size()) {
					return uniprotEvidenceArray;
				}
			}
		}

		return ret;
	}

	private ProteinEvidence getEvidence(ProteinGroup group, String acc) {
		if (group != null) {
			for (GroupableProtein groupableProtein : group) {
				if (groupableProtein.getAccession().equals(acc)) {
					return groupableProtein.getEvidence();
				}
			}
		}
		return null;
	}

	private void writeHeaders(FileWriter fw, CellType celltype) throws IOException {

		// write the header
		fw.write("NUM\t" + "VALID\t" + "ACC\t" + "Gene\t" + "protein description\t" + "PROTEIN_EVIDENCE\t" + "SCORE\t"
				+ "Known NE\t" + "Transmembrane region\t" + "Total SPC in " + celltype + "\t");

		List<Experiment> experimentList = getAllExperiments();

		if (celltype != null) {
			for (Experiment experiment : experimentList) {
				if (experiment.getCellType() == celltype) {
					// SPC replicates
					for (CellCompartment cellCompartment : CellCompartment.values()) {
						List<Replicate> replicates = experiment.getSortedReplicates();
						for (Replicate replicate : replicates) {
							Fractionation fractionation = replicate.getFractionation(cellCompartment);
							if (fractionation != null) {
								fw.write(fractionation.getName() + "_SPC\t");
							}
						}
					}

				}
			}
			for (Experiment experiment : experimentList) {
				if (experiment.getCellType() == celltype) {
					// NSAF replicates
					for (CellCompartment cellCompartment : CellCompartment.values()) {
						List<Replicate> replicates = experiment.getSortedReplicates();
						for (Replicate replicate : replicates) {
							Fractionation fractionation = replicate.getFractionation(cellCompartment);
							if (fractionation != null) {
								fw.write(fractionation.getName() + "_NSAF\t");
							}
						}
					}
				}
			}
		} else {
			throw new IllegalArgumentException("CellType null not supported.");
		}
		fw.write("\n");

	}

	public Set<String> getAllAccs(CellType cellType) throws IOException {
		if (proteinAccs == null) {
			proteinAccs = new HashSet<String>();
			List<Experiment> allExperiments = getAllExperiments();
			for (Experiment experiment : allExperiments) {
				if (experiment.getCellType() == cellType) {
					proteinAccs.addAll(experiment.getProteinAccs());
				}
			}
		}
		return proteinAccs;
	}

	public List<String> getAccs(String rawAcc) {
		List<String> accs = new ArrayList<String>();

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
				if (!accs.contains(acc)) {
					accs.add(acc);
				}
			}
		} else {
			if (!accs.contains(rawAcc)) {
				accs.add(rawAcc);
			}
		}
		return accs;

	}

	public String getProteinNameString(String rawAcc, CellType cellType) throws IOException {

		String proteinName = "";
		List<Boolean> validArray = filterAccessionsByEvidence(rawAcc, null, cellType);
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

	public String getGeneNameString(String rawAcc, CellType cellType) throws IOException {

		List<String> geneNames = new ArrayList<String>();
		List<Boolean> validArray = filterAccessionsByEvidence(rawAcc, null, cellType);
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

	public List<Experiment> getAllExperiments() {
		List<Experiment> list = new ArrayList<Experiment>();
		list.addAll(experimentsA);
		list.addAll(experimentsM);
		list.addAll(experimentsU);
		return list;
	}

	private ProteinGroup getGroup(String groupKey, CellType cellType) throws IOException {
		if (!groupsByRawAcc.containsKey(groupKey)) {
			List<ProteinGroup> proteinGroupsTotal = getProteinGroups(cellType);
			List<ProteinGroup> proteinGroups = proteinGroupsTotal.stream().filter(
					group -> group.getKey().contains(groupKey) && !group.getKey().toLowerCase().contains("reverse"))
					.collect(Collectors.toList());
			if (proteinGroups.size() > 1) {
				log.warn("More than one group  with key " + groupKey);
			}
			if (proteinGroups.isEmpty()) {
				proteinGroups = proteinGroupsTotal.stream().filter(group -> group.getKey().contains(groupKey))
						.collect(Collectors.toList());
				if (!proteinGroups.isEmpty()) {
					groupsByRawAcc.put(groupKey, proteinGroups.get(0));
				}
			}
			if (!proteinGroups.isEmpty()) {
				groupsByRawAcc.put(groupKey, proteinGroups.get(0));
			} else {
				log.warn("Group not found for " + groupKey + " in celltype:" + cellType);
			}
		}
		return groupsByRawAcc.get(groupKey);
	}

	private String getFilteredProteinEvidences(String rawAcc, CellType cellType) throws IOException {
		if (rawAcc.contains("Q99M74")) {
			log.warn(rawAcc);
		}
		List<Boolean> validArray = filterAccessionsByEvidence(rawAcc, null, cellType);
		List<ProteinEvidence> evidences = new ArrayList<ProteinEvidence>();
		int index = 0;
		ProteinGroup group = getGroup(rawAcc, null);
		for (String acc : getAccs(rawAcc)) {
			if (!validArray.get(index++)) {
				continue;
			}
			ProteinEvidence evidence = getEvidence(group, acc);
			if (!evidences.contains(evidence)) {
				evidences.add(evidence);
			}
		}
		StringBuilder sb = new StringBuilder();

		for (ProteinEvidence evidence : evidences) {
			if (!"".equals(sb.toString())) {
				sb.append(Constants.SEPARATOR);
			}
			if (evidence != null) {
				sb.append(evidence);
			} else {

				for (GroupableProtein groupableProtein : group) {
					log.warn(groupableProtein.getAccession());
				}
				sb.append("");
			}
		}
		return sb.toString();
	}

	private String getProteinEvidences(String rawAcc) throws IOException {

		List<String> proteinAccs = getAccs(rawAcc);

		ProteinGroup group = getGroup(rawAcc, null);
		StringBuilder sb = new StringBuilder();
		if (group != null) {
			for (String proteinAcc : proteinAccs) {
				ProteinEvidence evidence = getEvidence(group, proteinAcc);
				if (!"".equals(sb.toString())) {
					sb.append(Constants.SEPARATOR);
				}
				if (evidence != null) {
					sb.append(evidence);
				} else {
					log.warn(" I didnt find " + proteinAcc + "  in group  " + group);
					for (GroupableProtein groupableProtein : group) {
						log.warn(groupableProtein.getAccession());
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

	private Filter getKeratinFilter() {

		if (Constants.KERATIN_FILTER && keratinFilter == null) {
			keratinFilter = new KeratinFilter("ASDF");
		}
		return keratinFilter;
	}

	public boolean isValid(String rawAcc, int totalSPC) throws IOException {
		if (totalSPC < Constants.MIN_TOTAL_SPC) {
			return false;
		}
		final Filter goFilter2 = getGOFilter();
		final Filter keratinFilter = getKeratinFilter();
		if (goFilter2 != null || keratinFilter != null) {
			final ProteinGroup group = getGroup(rawAcc, null);
			if (group != null) {
				for (GroupableProtein groupableProtein : group) {
					if (goFilter2 != null && !goFilter2.isValid((Protein) groupableProtein)) {
						return false;
					}
					if (keratinFilter != null && !keratinFilter.isValid((Protein) groupableProtein)) {
						return false;
					}
				}
				return true;
			}
		}
		return true;
	}

}
