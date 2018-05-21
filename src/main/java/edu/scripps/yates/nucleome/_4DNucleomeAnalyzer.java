package edu.scripps.yates.nucleome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
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
import edu.scripps.yates.nucleome.model.Wash;
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
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.hash.THashSet;

public class _4DNucleomeAnalyzer {
	private final static Logger log = Logger.getLogger(_4DNucleomeAnalyzer.class);

	public static void main(String[] args) {
		_4DNucleomeAnalyzer analyzer;
		try {

			analyzer = new _4DNucleomeAnalyzer(args[0]);

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
			Constants.writeCoverageFile = true;
			UniprotProteinRetriever.enableCache = true;
			scoringFunction = new ScoringFunctionByNE_NSAF_Ratios(analyzer);
			// scoringFunction = new ScoringFunctionByNE_NSAF_Points(analyzer);
			////////////////////////////////////////////////////////////

			analyzer.run();

			System.out.println("DONE");
			System.exit(0);
		} catch (final IOException e) {
			e.printStackTrace();
			System.err.println("ERROR: " + e.getMessage());
		}
		System.exit(-1);
	}

	private final List<Experiment> experimentsU = new ArrayList<>();
	private final List<Experiment> experimentsA = new ArrayList<>();
	private final List<Experiment> experimentsM = new ArrayList<>();
	protected static ScoringFunction scoringFunction;
	private final String hostName = "garibaldi.scripps.edu";;
	private final String userName = "salvador";
	private final String pass;
	private final String remotefileName = "DTASelect-filter.txt";
	protected static String datasetsPathsFile = "z:\\share\\Salva\\data\\4D_Nucleome\\SwissProt_1FDR\\SwissProt_1FDR_data_paths.txt";

	private static String datasetsPhosphoPathsFile = "z:\\share\\Salva\\data\\4D_Nucleome\\datasets_paths_phospho.txt";
	protected static File outputFolder = new File(new File(datasetsPathsFile).getParent() + File.separator + "output");

	private final Map<CellType, List<Pair<String, Double>>> scoresByCellType = new THashMap<CellType, List<Pair<String, Double>>>();
	private GOFilter goFilter;
	private KeratinFilter keratinFilter;
	private int fileNum = 1;
	protected List<ProteinGroup> proteinGroups;

	protected Set<String> proteinAccs;
	protected TObjectIntHashMap<String> totalSPCs = new TObjectIntHashMap<String>();
	protected Map<String, ProteinGroup> groupsByRawAcc = new THashMap<String, ProteinGroup>();
	protected Set<GroupableProtein> groupableProteins = new THashSet<GroupableProtein>();

	public _4DNucleomeAnalyzer(String pass) throws IOException {
		this.pass = pass;
	}

	/**
	 * Returns a list of pairs (Protein ACC - Score), sorted by the score
	 *
	 * @return
	 * @throws IOException
	 */

	private List<Pair<String, Double>> calculateScoresFromGroups(CellType celltype, Wash wash) throws IOException {
		log.info("Calculating scores for " + getAllAccs(celltype).size() + " proteins in " + celltype + " Wash: "
				+ wash);
		final List<Pair<String, Double>> scores = new ArrayList<Pair<String, Double>>();
		final List<ProteinGroup> proteinGroups = getProteinGroups(celltype, wash);

		for (final ProteinGroup proteinGroup : proteinGroups) {
			if (proteinGroup.getEvidence() == ProteinEvidence.NONCONCLUSIVE) {
				continue;
			}
			final String rawAcc = proteinGroup.getKey();

			final List<String> filteredAcessions = new ArrayList<String>();
			final String filteredAcc = getAccessionStringByEvidence(rawAcc, proteinGroup, celltype, wash);
			if (filteredAcc.contains(",")) {
				final String[] split = filteredAcc.split(",");
				for (final String acc : split) {
					filteredAcessions.add(acc);
				}
			} else {
				filteredAcessions.add(filteredAcc);
			}
			final double score = scoringFunction.getScore(filteredAcessions, celltype, wash);

			final Pair<String, Double> pair = new Pair<>(filteredAcc, score);
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

					final double score1 = o1.getSecondElement();
					final int totalSPC1 = getTotalSPC(rawAcc1, celltype, wash);
					final boolean valid1 = isValid(rawAcc1, totalSPC1) ? true : false;
					//
					String rawAcc2 = o2.getFirstelement();
					if (rawAcc2.contains("[")) {
						rawAcc2 = rawAcc2.substring(0, rawAcc2.indexOf("["));
					}
					final double score2 = o2.getSecondElement();
					final int totalSPC2 = getTotalSPC(rawAcc2, celltype, wash);
					final boolean valid2 = isValid(rawAcc2, totalSPC2) ? true : false;
					//

					if (valid1 && !valid2) {
						return -1;
					} else if (!valid1 && valid2) {
						return 1;
					} else {
						final int scoreComparison = Double.compare(score2, score1);
						if (scoreComparison == 0) {
							// by total SPC
							return Integer.compare(totalSPC2, totalSPC1);
						} else {
							return scoreComparison;
						}
					}

				} catch (final IOException e) {
					return 0;
				}
			}
		});
		log.info(scores.size() + " scores calculated");

		scoresByCellType.put(celltype, scores);
		return scores;
	}

	private List<ProteinGroup> getProteinGroups(CellType cellType, Wash wash) throws IOException {
		if (proteinGroups == null) {
			final PAnalyzer panalyzer = new PAnalyzer(true);
			final Set<GroupableProtein> groupableProteins = getAllGroupableProteins(cellType, wash);

			log.info("Grouping " + groupableProteins.size() + " proteins from " + cellType + " ans wash " + wash);
			final List<ProteinGroup> proteinGroupsTMP = panalyzer.run(groupableProteins);
			log.info(proteinGroupsTMP.size() + " protein groups");
			proteinGroups = new ArrayList<>();
			for (final ProteinGroup proteinGroup : proteinGroupsTMP) {
				if (proteinGroup.getEvidence() != ProteinEvidence.NONCONCLUSIVE) {
					proteinGroups.add(proteinGroup);
				}
			}
			log.info(proteinGroups.size() + " protein groups after removing NonConclusive proteins");
		}
		return proteinGroups;
	}

	private Set<GroupableProtein> getAllGroupableProteins(CellType cellType, Wash wash) throws IOException {
		if (groupableProteins.isEmpty()) {
			for (final Experiment experiment : getAllExperiments()) {
				if (cellType != null && !experiment.getCellType().equals(cellType)) {
					continue;
				}
				if (wash != null && !experiment.getWash().equals(wash)) {
					continue;
				}
				for (final Protein protein : experiment.getProteins()) {
					groupableProteins.add(protein);
				}
			}

		}
		return groupableProteins;
	}

	private void loadDatasets() throws IOException {
		log.info("Loading datasets");
		final long t1 = System.currentTimeMillis();

		experimentsU.clear();
		experimentsA.clear();
		experimentsM.clear();

		final DataPaths dataPaths = new DataPaths(Constants.DATASET_PATHS_FILE);
		// U (N, Ne, C)

		// FDR 1%

		final String[] u2Files = dataPaths.getFiles("U2");
		final String[] u3Files = dataPaths.getFiles("U3");
		final String[] u4Files = dataPaths.getFiles("U41");
		final String[] u42Files = dataPaths.getFiles("U42");
		// String[] u5Files = dataPaths.getFiles("U5");
		final String[] a2Files = dataPaths.getFiles("A2");
		final String[] a4Files = dataPaths.getFiles("A41");
		final String[] a42Files = dataPaths.getFiles("A42");
		final String[] m1Files = dataPaths.getFiles("M11");
		final String[] m12Files = dataPaths.getFiles("M12");
		final String[] m3Files = dataPaths.getFiles("M31");
		final String[] m32Files = dataPaths.getFiles("M32");

		final Experiment experimentU2 = new Experiment("U2", CellType.U);
		experimentsU.add(experimentU2);
		experimentU2.addReplicate(1, CellType.U, CellCompartment.N, getRemoteFile(u2Files[0]));
		experimentU2.addReplicate(1, CellType.U, CellCompartment.NE, getRemoteFile(u2Files[1]));
		experimentU2.addReplicate(1, CellType.U, CellCompartment.CM, getRemoteFile(u2Files[2]));

		final Experiment experimentU3 = new Experiment("U3", CellType.U);
		experimentsU.add(experimentU3);
		experimentU3.addReplicate(1, CellType.U, CellCompartment.N, getRemoteFile(u3Files[0]));
		experimentU3.addReplicate(1, CellType.U, CellCompartment.NE, getRemoteFile(u3Files[1]));
		experimentU3.addReplicate(1, CellType.U, CellCompartment.CM, getRemoteFile(u3Files[2]));

		final Experiment experimentU4 = new Experiment("U4", CellType.U);
		experimentsU.add(experimentU4);
		experimentU4.addReplicate(1, CellType.U, CellCompartment.N, getRemoteFile(u4Files[0]));
		experimentU4.addReplicate(1, CellType.U, CellCompartment.NE, getRemoteFile(u4Files[1]));
		experimentU4.addReplicate(1, CellType.U, CellCompartment.CM, getRemoteFile(u4Files[2]));
		experimentU4.addReplicate(2, CellType.U, CellCompartment.N, getRemoteFile(u42Files[0]));
		experimentU4.addReplicate(2, CellType.U, CellCompartment.NE, getRemoteFile(u42Files[1]));
		experimentU4.addReplicate(2, CellType.U, CellCompartment.CM, getRemoteFile(u42Files[2]));

		// Experiment experimentU5 = new Experiment("U5", CellType.U);
		// experimentsU.add(experimentU5);
		// experimentU5.addReplicate(1, CellType.U, CellCompartment.NE,
		// getRemoteFile(u5Files[0]));
		// experimentU5.addReplicate(1, CellType.U, CellCompartment.C,
		// getRemoteFile(u5Files[1]));

		if (!Constants.TESTING) {
			final Experiment experimentA3 = new Experiment("A2", CellType.A);
			experimentsA.add(experimentA3);
			experimentA3.addReplicate(1, CellType.A, CellCompartment.N, getRemoteFile(a2Files[0]));
			experimentA3.addReplicate(1, CellType.A, CellCompartment.NE, getRemoteFile(a2Files[1]));
			experimentA3.addReplicate(1, CellType.A, CellCompartment.CM, getRemoteFile(a2Files[2]));

			final Experiment experimentA4 = new Experiment("A4", CellType.A);
			experimentsA.add(experimentA4);
			experimentA4.addReplicate(1, CellType.A, CellCompartment.N, getRemoteFile(a4Files[0]));
			experimentA4.addReplicate(1, CellType.A, CellCompartment.NE, getRemoteFile(a4Files[1]));
			experimentA4.addReplicate(1, CellType.A, CellCompartment.CM, getRemoteFile(a4Files[2]));
			experimentA4.addReplicate(2, CellType.A, CellCompartment.N, getRemoteFile(a42Files[0]));
			experimentA4.addReplicate(2, CellType.A, CellCompartment.NE, getRemoteFile(a42Files[1]));
			experimentA4.addReplicate(2, CellType.A, CellCompartment.CM, getRemoteFile(a42Files[2]));

			final Experiment experimentM1 = new Experiment("M1", CellType.M);
			experimentsM.add(experimentM1);
			experimentM1.addReplicate(1, CellType.M, CellCompartment.N, getRemoteFile(m1Files[0]));
			experimentM1.addReplicate(1, CellType.M, CellCompartment.NE, getRemoteFile(m1Files[1]));
			experimentM1.addReplicate(1, CellType.M, CellCompartment.CM, getRemoteFile(m1Files[2]));
			experimentM1.addReplicate(2, CellType.M, CellCompartment.N, getRemoteFile(m12Files[0]));
			experimentM1.addReplicate(2, CellType.M, CellCompartment.NE, getRemoteFile(m12Files[1]));
			experimentM1.addReplicate(2, CellType.M, CellCompartment.CM, getRemoteFile(m12Files[2]));

			final Experiment experimentM3 = new Experiment("M3", CellType.M);
			experimentsM.add(experimentM3);

			experimentM3.addReplicate(1, CellType.M, CellCompartment.N, getRemoteFile(m3Files[0]));
			experimentM3.addReplicate(1, CellType.M, CellCompartment.NE, getRemoteFile(m3Files[1]));
			experimentM3.addReplicate(1, CellType.M, CellCompartment.CM, getRemoteFile(m3Files[2]));
			experimentM3.addReplicate(2, CellType.M, CellCompartment.N, getRemoteFile(m32Files[0]));
			experimentM3.addReplicate(2, CellType.M, CellCompartment.NE, getRemoteFile(m32Files[1]));
			experimentM3.addReplicate(2, CellType.M, CellCompartment.CM, getRemoteFile(m32Files[2]));

		}

		final long t2 = System.currentTimeMillis();
		log.info("It took " + DatesUtil.getDescriptiveTimeFromMillisecs(t2 - t1));

	}

	private void loadDatasetsXi() throws IOException {
		log.info("Loading datasets");
		final long t1 = System.currentTimeMillis();

		experimentsU.clear();
		experimentsA.clear();
		experimentsM.clear();

		final DataPaths dataPaths = new DataPaths(Constants.DATASET_PATHS_FILE);
		// U (N, Ne, C)

		// FDR 1%

		final String[] u2Files = dataPaths.getFiles("U2");
		final String[] u3Files = dataPaths.getFiles("U3");
		final String[] u4Files = dataPaths.getFiles("U41");
		final String[] u42Files = dataPaths.getFiles("U42");
		// String[] u5Files = dataPaths.getFiles("U5");
		final String[] a2Files = dataPaths.getFiles("A2");
		final String[] a4Files = dataPaths.getFiles("A41");
		final String[] a42Files = dataPaths.getFiles("A42");
		final String[] m1Files = dataPaths.getFiles("M11");
		final String[] m12Files = dataPaths.getFiles("M12");
		final String[] m3Files = dataPaths.getFiles("M31");
		final String[] m32Files = dataPaths.getFiles("M32");

		final Experiment experimentU2 = new Experiment("U2", CellType.U);
		experimentsU.add(experimentU2);
		experimentU2.addReplicate(1, CellType.U, CellCompartment.N, getRemoteFile(u2Files[0]));
		experimentU2.addReplicate(1, CellType.U, CellCompartment.NE, getRemoteFile(u2Files[1]));
		experimentU2.addReplicate(1, CellType.U, CellCompartment.CM, getRemoteFile(u2Files[2]));

		final Experiment experimentU3 = new Experiment("U3", CellType.U);
		experimentsU.add(experimentU3);
		experimentU3.addReplicate(1, CellType.U, CellCompartment.N, getRemoteFile(u3Files[0]));
		experimentU3.addReplicate(1, CellType.U, CellCompartment.NE, getRemoteFile(u3Files[1]));
		experimentU3.addReplicate(1, CellType.U, CellCompartment.CM, getRemoteFile(u3Files[2]));

		final Experiment experimentU4 = new Experiment("U4", CellType.U);
		experimentsU.add(experimentU4);
		experimentU4.addReplicate(1, CellType.U, CellCompartment.N, getRemoteFile(u4Files[0]));
		experimentU4.addReplicate(1, CellType.U, CellCompartment.NE, getRemoteFile(u4Files[1]));
		experimentU4.addReplicate(1, CellType.U, CellCompartment.CM, getRemoteFile(u4Files[2]));
		experimentU4.addReplicate(2, CellType.U, CellCompartment.N, getRemoteFile(u42Files[0]));
		experimentU4.addReplicate(2, CellType.U, CellCompartment.NE, getRemoteFile(u42Files[1]));
		experimentU4.addReplicate(2, CellType.U, CellCompartment.CM, getRemoteFile(u42Files[2]));

		// Experiment experimentU5 = new Experiment("U5", CellType.U);
		// experimentsU.add(experimentU5);
		// experimentU5.addReplicate(1, CellType.U, CellCompartment.NE,
		// getRemoteFile(u5Files[0]));
		// experimentU5.addReplicate(1, CellType.U, CellCompartment.C,
		// getRemoteFile(u5Files[1]));

		if (!Constants.TESTING) {
			final Experiment experimentA3 = new Experiment("A2", CellType.A);
			experimentsA.add(experimentA3);
			experimentA3.addReplicate(1, CellType.A, CellCompartment.N, getRemoteFile(a2Files[0]));
			experimentA3.addReplicate(1, CellType.A, CellCompartment.NE, getRemoteFile(a2Files[1]));
			experimentA3.addReplicate(1, CellType.A, CellCompartment.CM, getRemoteFile(a2Files[2]));

			final Experiment experimentA4 = new Experiment("A4", CellType.A);
			experimentsA.add(experimentA4);
			experimentA4.addReplicate(1, CellType.A, CellCompartment.N, getRemoteFile(a4Files[0]));
			experimentA4.addReplicate(1, CellType.A, CellCompartment.NE, getRemoteFile(a4Files[1]));
			experimentA4.addReplicate(1, CellType.A, CellCompartment.CM, getRemoteFile(a4Files[2]));
			experimentA4.addReplicate(2, CellType.A, CellCompartment.N, getRemoteFile(a42Files[0]));
			experimentA4.addReplicate(2, CellType.A, CellCompartment.NE, getRemoteFile(a42Files[1]));
			experimentA4.addReplicate(2, CellType.A, CellCompartment.CM, getRemoteFile(a42Files[2]));

			final Experiment experimentM1 = new Experiment("M1", CellType.M);
			experimentsM.add(experimentM1);
			experimentM1.addReplicate(1, CellType.M, CellCompartment.N, getRemoteFile(m1Files[0]));
			experimentM1.addReplicate(1, CellType.M, CellCompartment.NE, getRemoteFile(m1Files[1]));
			experimentM1.addReplicate(1, CellType.M, CellCompartment.CM, getRemoteFile(m1Files[2]));
			experimentM1.addReplicate(2, CellType.M, CellCompartment.N, getRemoteFile(m12Files[0]));
			experimentM1.addReplicate(2, CellType.M, CellCompartment.NE, getRemoteFile(m12Files[1]));
			experimentM1.addReplicate(2, CellType.M, CellCompartment.CM, getRemoteFile(m12Files[2]));

			final Experiment experimentM3 = new Experiment("M3", CellType.M);
			experimentsM.add(experimentM3);

			experimentM3.addReplicate(1, CellType.M, CellCompartment.N, getRemoteFile(m3Files[0]));
			experimentM3.addReplicate(1, CellType.M, CellCompartment.NE, getRemoteFile(m3Files[1]));
			experimentM3.addReplicate(1, CellType.M, CellCompartment.CM, getRemoteFile(m3Files[2]));
			experimentM3.addReplicate(2, CellType.M, CellCompartment.N, getRemoteFile(m32Files[0]));
			experimentM3.addReplicate(2, CellType.M, CellCompartment.NE, getRemoteFile(m32Files[1]));
			experimentM3.addReplicate(2, CellType.M, CellCompartment.CM, getRemoteFile(m32Files[2]));

		}

		final long t2 = System.currentTimeMillis();
		log.info("It took " + DatesUtil.getDescriptiveTimeFromMillisecs(t2 - t1));

	}

	protected File getRemoteFile(Pair<String, String> remotePath) throws IOException {
		final File localFile = new File(outputFolder + File.separator + "data_files" + File.separator
				+ FilenameUtils.getBaseName(remotePath.getSecondElement() + remotePath.getFirstelement()) + "."
				+ FilenameUtils.getExtension(remotefileName));
		if (localFile.exists() && localFile.length() > 0) {
			return localFile;
		}
		final RemoteSSHFileReference ret = new RemoteSSHFileReference(hostName, userName, pass, remotefileName, null);
		ret.setRemotePath(remotePath.getSecondElement());

		ret.setOutputFile(localFile);
		final File remoteFile = ret.getRemoteFile();
		return remoteFile;
	}

	protected File getRemoteFile(String remotePath) throws IOException {
		final File localFile = new File(outputFolder + File.separator + "data_files" + File.separator
				+ FilenameUtils.getBaseName(remotePath + remotefileName) + "_" + fileNum++ + "."
				+ FilenameUtils.getExtension(remotefileName));
		if (localFile.exists() && localFile.length() > 0) {
			return localFile;
		}
		final RemoteSSHFileReference ret = new RemoteSSHFileReference(hostName, userName, pass, remotefileName, null);
		ret.setRemotePath(remotePath);

		ret.setOutputFile(localFile);
		final File remoteFile = ret.getRemoteFile();
		return remoteFile;
	}

	public void run() throws IOException {
		// choose scoring function
		final long t1 = System.currentTimeMillis();
		try {
			// load data
			loadDatasets();
			// print scores for each cell type
			final CellType[] values = CellType.values();
			for (final CellType cellType : values) {
				for (final Wash wash : Wash.values()) {
					// restart field variables
					proteinAccs = null;
					proteinGroups = null;
					totalSPCs.clear();
					groupsByRawAcc.clear();
					groupableProteins.clear();

					// annotate proteins with uniprot
					annotateProteins(cellType);
					if (!Constants.writeCoverageFile) {
						writeScoreDistributions(cellType, wash);
					}

					// writeScoreDistributions(cellType, DataType.NSAF, true &&
					// Constants.printRatios);
					// writeScoreDistributions(cellType, DataType.PEPC, false &&
					// Constants.printRatios);
				}
			}
			if (Constants.writeCoverageFile) {
				writeCoverageFile();
			}
			// print scores for all celltypes together
			// writeScoreDistributions(null, DataType.NSAF, true &&
			// Constants.printRatios);
			if (Constants.writeCombinedDistribution) {
				writeScoreDistributions(null, null);
			}
			// writeScoreDistributions(null, DataType.PEPC, false &&
			// Constants.printRatios);

			if (Constants.compareScores) {
				// compare the scores between U and A
				final PairComparisonReport comparisonReportUA = compareScores(CellType.U, null, CellType.A, null);
				comparisonReportUA.printToFile(
						new File(outputFolder.getAbsolutePath() + File.separator + "U_vs_A_comparison.txt"));
				// compare the scores between U and M
				final PairComparisonReport comparisonReportUM = compareScores(CellType.U, null, CellType.M, null);
				comparisonReportUM.printToFile(
						new File(outputFolder.getAbsolutePath() + File.separator + "U_vs_M_comparison.txt"));
				// compare the scores between A and M
				final PairComparisonReport comparisonReportAM = compareScores(CellType.A, null, CellType.M, null);
				comparisonReportAM.printToFile(
						new File(outputFolder.getAbsolutePath() + File.separator + "A_vs_M_comparison.txt"));
				// compare the scores between U and A and M
				final TripleComparisonReport comparisonReportUAM = compareScores(CellType.U, null, CellType.A, null,
						CellType.M, null);
				comparisonReportUAM.printToFile(
						new File(outputFolder.getAbsolutePath() + File.separator + "U_vs_A_vs_M_comparison.txt"));
			}
		} finally {
			log.info("It took " + DatesUtil.getDescriptiveTimeFromMillisecs(System.currentTimeMillis() - t1));
		}
	}

	private void writeCoverageFile() throws IOException {
		final List<ProteinGroup> proteinGroups = getProteinGroups(null, null);
		final FileWriter fw = new FileWriter(outputFolder + File.separator + "coverages.txt");
		// header
		fw.write("\t");
		final List<Experiment> allExperiments = getAllExperiments();
		for (final Experiment experiment : allExperiments) {
			final List<Replicate> replicates = experiment.getReplicates();
			for (final Replicate replicate : replicates) {
				for (final CellCompartment cellCompartment : CellCompartment.values()) {
					final Fractionation fractionation = replicate.getFractionation(cellCompartment);
					if (fractionation != null) {
						fw.write(fractionation.getName() + "_%cov\t");
					}
				}
			}
		}
		fw.write("\n");
		for (final ProteinGroup proteinGroup : proteinGroups) {
			if (proteinGroup.getEvidence() == ProteinEvidence.NONCONCLUSIVE) {
				continue;
			}
			final String rawAcc = proteinGroup.getKey();

			final List<String> filteredAcessions = new ArrayList<String>();
			final String filteredAcc = getAccessionStringByEvidence(rawAcc, proteinGroup, null, null);
			if (filteredAcc.contains(",")) {
				final String[] split = filteredAcc.split(",");
				for (final String acc : split) {
					filteredAcessions.add(acc);
				}
			} else {
				filteredAcessions.add(filteredAcc);
			}
			fw.write(proteinGroup.getKey() + "\t");

			// look into all the experiments

			for (final Experiment experiment : allExperiments) {
				final List<Replicate> replicates = experiment.getReplicates();
				for (final Replicate replicate : replicates) {
					for (final CellCompartment cellCompartment : CellCompartment.values()) {
						final Fractionation fractionation = replicate.getFractionation(cellCompartment);
						if (fractionation != null) {
							final String coverage = fractionation.getCoverage(filteredAcessions, true);
							fw.write(coverage + "\t");
						}
					}
				}

			}

			fw.write("\n");

		}

		fw.close();

	}

	protected void annotateProteins(CellType cellType) throws IOException {

		final Set<String> uniprotAccs = new THashSet<String>();
		for (final String rawAcc : getAllAccs(cellType)) {
			final String uniprotAcc = FastaParser.getUniProtACC(rawAcc);
			if (uniprotAcc == null) {
				log.warn("Not uniprot acc: " + rawAcc);
				continue;
			} else {
				uniprotAccs.add(uniprotAcc);
			}
		}
		final UniprotProteinRetriever upr = Constants.upr;
		final Map<String, Protein> annotatedProteins = upr.getAnnotatedProteins(uniprotAccs);
		log.info(annotatedProteins.size() + " proteins annotated out of " + uniprotAccs.size());
	}

	private PairComparisonReport compareScores(CellType cellType1, Wash wash1, CellType cellType2, Wash wash2) {
		log.info(" Comparing " + cellType1 + " with " + cellType2);
		final List<Pair<String, Double>> scores1 = scoresByCellType.get(cellType1);
		final Map<String, Double> scoreMap1 = getScoreMap(scores1);
		final List<Pair<String, Double>> scores2 = scoresByCellType.get(cellType2);
		final Map<String, Double> scoreMap2 = getScoreMap(scores2);
		return new PairComparisonReport(this, cellType1, wash1, scoreMap1, cellType2, wash2, scoreMap2);
	}

	private TripleComparisonReport compareScores(CellType cellType1, Wash wash1, CellType cellType2, Wash wash2,
			CellType cellType3, Wash wash3) {
		final List<Pair<String, Double>> scores1 = scoresByCellType.get(cellType1);
		final Map<String, Double> scoreMap1 = getScoreMap(scores1);
		final List<Pair<String, Double>> scores2 = scoresByCellType.get(cellType2);
		final Map<String, Double> scoreMap2 = getScoreMap(scores2);
		final List<Pair<String, Double>> scores3 = scoresByCellType.get(cellType3);
		final Map<String, Double> scoreMap3 = getScoreMap(scores3);
		return new TripleComparisonReport(this, cellType1, wash1, scoreMap1, cellType2, wash2, scoreMap2, cellType3,
				wash3, scoreMap3);
	}

	private Map<String, Double> getScoreMap(List<Pair<String, Double>> scoreList) {

		final Map<String, Double> ret = new THashMap<String, Double>();
		for (final Pair<String, Double> pair : scoreList) {
			ret.put(pair.getFirstelement(), pair.getSecondElement());
		}
		return ret;
	}

	protected void writeScoreDistributions(CellType celltype, Wash wash) throws IOException {
		if (!Constants.printScoreDistributions) {
			return;
		}
		final List<Pair<String, Double>> scores = calculateScoresFromGroups(celltype, wash);

		log.info("Printing scores for " + celltype + " to file");
		String cellTypeName = "UAM";
		if (celltype != null) {
			cellTypeName = celltype.name();
		}
		final String formatedDate = DateFormatUtils.format(new Date(), "yyyy-MM-dd_HH-mm");
		final String pathname = outputFolder.getAbsolutePath() + File.separator + formatedDate + "_" + cellTypeName
				+ "_" + wash + "_scores_distribution_" + scoringFunction.getName() + ".txt";

		final File scoreFileOutput = new File(pathname);
		if (!scoreFileOutput.getParentFile().exists()) {
			log.info("Creating file " + scoreFileOutput.getAbsolutePath());
			scoreFileOutput.getParentFile().mkdirs();
		}
		if (scoreFileOutput.exists()) {
			log.info("Overriding file " + scoreFileOutput.getAbsolutePath());
		}
		FileWriter fw = null;
		try {
			final ProgressCounter counter = new ProgressCounter(scores.size(), ProgressPrintingType.PERCENTAGE_STEPS,
					0);
			fw = new FileWriter(scoreFileOutput);
			writeHeaders(fw, celltype, wash);
			int num = 1;
			for (final Pair<String, Double> pair : scores) {
				counter.increment();
				final String percentage = counter.printIfNecessary();
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
					final String rawAccString = getAccessionStringByEvidence(rawAcc, null, celltype, wash);
					if (rawAccString.equals("Q91W53")) {
						log.info(rawAcc);
					}
					final List<String> filteredAccessions = getAccs(rawAccString);
					String geneNameString = null;
					// TODO
					// if (false) {
					geneNameString = getGeneNameString(rawAcc, celltype, wash);
					// }

					final String proteinNameString = getProteinNameString(rawAcc, celltype, wash);
					final double score = pair.getSecondElement();
					// end testing
					final int totalSPC = getTotalSPC(rawAcc, celltype, wash);
					final String valid = isValid(rawAcc, totalSPC) ? "VALID" : "FILTERED";

					fw.write(num++ + "\t" + valid + "\t" + rawAccString + "\t" + geneNameString + "\t"
							+ proteinNameString + "\t" + getFilteredProteinEvidences(rawAcc, celltype, wash) + "\t"
							+ score + "\t" + ControlNE.isControl(rawAcc) + "\t"
							+ getTransmembraneRegion(rawAcc, celltype, wash) + "\t" + totalSPC + "\t");
					if (celltype != null) {
						for (final Experiment experiment : getAllExperiments()) {
							if (experiment.getCellType() == celltype) {
								if (wash == null || wash == experiment.getWash()) {
									// print the averages
									for (final CellCompartment cellCompartment : CellCompartment.values()) {
										final List<Replicate> replicates = experiment.getSortedReplicates();
										// SPCs
										for (final Replicate replicate : replicates) {
											final Fractionation fractionation = replicate
													.getFractionation(cellCompartment);
											if (fractionation != null) {
												final int repSPC = fractionation.getSpectralCount(filteredAccessions,
														true);
												fw.write(repSPC + "\t");
											}
										}
									}
								}
							}
						}
						for (final Experiment experiment : getAllExperiments()) {
							if (experiment.getCellType() == celltype) {
								if (wash == null || wash == experiment.getWash()) {
									// print the averages
									for (final CellCompartment cellCompartment : CellCompartment.values()) {
										final List<Replicate> replicates = experiment.getSortedReplicates();
										// NSAFs
										for (final Replicate replicate : replicates) {
											final Fractionation fractionation = replicate
													.getFractionation(cellCompartment);
											if (fractionation != null) {
												final double repNSAF = fractionation.getAverageNSAF(filteredAccessions,
														true);
												fw.write(repNSAF + "\t");
											}
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
		} catch (final Exception e) {
			e.printStackTrace();
			log.error(e.getMessage());
		} finally {
			if (fw != null) {
				fw.close();
			}
			log.info("File created at " + scoreFileOutput.getAbsolutePath());
		}

	}

	public int getTotalSPC(String rawAcc, CellType cellType, Wash wash) throws IOException {
		if (!totalSPCs.containsKey(rawAcc)) {

			final Set<String> individualAccs = new THashSet<String>();
			if (rawAcc.contains(",")) {
				// because if this has a comma it is because
				// it is an indistinguisable group
				final String[] split = rawAcc.split(",");
				for (final String string : split) {
					individualAccs.add(string);
				}
			} else {
				individualAccs.add(rawAcc);
			}

			final Set<GroupablePSM> psms = new THashSet<GroupablePSM>();

			final Set<GroupableProtein> proteins = getAllGroupableProteins(cellType, wash);
			for (final GroupableProtein groupableProtein : proteins) {
				if (rawAcc.contains(groupableProtein.getAccession())) {
					psms.addAll(groupableProtein.getGroupablePSMs());
				}
			}
			final int numPSMs = psms.size();
			int numPSMsTMP = 0;
			for (final Experiment experiment : getAllExperiments()) {
				if (experiment.getCellType() == cellType) {
					if (wash == null || experiment.getWash() == wash) {
						// print the averages
						for (final CellCompartment cellCompartment : CellCompartment.values()) {
							final List<Replicate> replicates = experiment.getSortedReplicates();
							// SPCs
							for (final Replicate replicate : replicates) {
								final Fractionation fractionation = replicate.getFractionation(cellCompartment);
								if (fractionation != null) {
									for (final String acc : individualAccs) {

										final int spectralCount = fractionation.getSpectralCount(acc, true);
										numPSMsTMP += spectralCount;
										// log.info(acc + " " + spectralCount +
										// " in
										// " + fractionation.getName());
										if (spectralCount > 0) {
											// we check all of them but we only
											// count one
											// in some dtaselects may appear the
											// first one and in other may appear
											// the
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
			}
			if (numPSMs != numPSMsTMP) {
				log.warn("cuidado ");
			}
			totalSPCs.put(rawAcc, numPSMs);
		}
		return totalSPCs.get(rawAcc);
	}

	private boolean getTransmembraneRegion(String rawAcc, CellType cellType, Wash wash) throws IOException {

		final List<Boolean> validArray = filterAccessionsByEvidence(rawAcc, null, cellType, wash);
		int index = 0;
		for (final String acc : getAccs(rawAcc)) {
			if (!validArray.get(index++)) {
				continue;
			}
			final Map<String, Protein> annotatedProtein = Constants.upr.getAnnotatedProtein(acc);
			if (annotatedProtein.containsKey(acc)) {
				final Protein protein2 = annotatedProtein.get(acc);
				if (protein2 != null) {
					final Set<ProteinAnnotation> annotations = protein2.getAnnotations();
					for (final ProteinAnnotation proteinAnnotation : annotations) {
						if (proteinAnnotation.getAnnotationType().getKey().equals("transmembrane region")) {
							return true;
						}
					}

				}
			}
		}
		return false;
	}

	private String getAccessionStringByEvidence(String rawAcc, ProteinGroup proteinGroup, CellType cellType, Wash wash)
			throws IOException {
		final List<Boolean> validArray = filterAccessionsByEvidence(rawAcc, proteinGroup, cellType, wash);
		int index = 0;
		final StringBuilder sb = new StringBuilder();
		for (final String acc : getAccs(rawAcc)) {
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

	private List<Boolean> filterAccessionsByEvidence(String rawAcc, ProteinGroup group, CellType cellType, Wash wash)
			throws IOException {

		final List<String> accs = getAccs(rawAcc);
		final Map<String, Entry> uniprotEntries = Constants.upr.getAnnotatedUniprotEntries(accs);
		if (group == null) {
			group = getGroup(rawAcc, cellType, wash);
		}
		final List<Boolean> ret = new ArrayList<>();
		final List<Boolean> groupEvidenceArray = new ArrayList<>();
		final List<Boolean> uniprotEvidenceArray = new ArrayList<>();
		// only swissprot is valid
		boolean thereIsASwissProt = false;
		boolean thereIsAConclusiveProt = false;
		for (final String acc : accs) {
			boolean valid = false;
			try {

				final ProteinEvidence evidence = getEvidence(group, acc);
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
							final String dataset = protein2.getDataset();
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
		for (final Boolean uniprotEvidence : uniprotEvidenceArray) {
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
			for (final GroupableProtein groupableProtein : group) {
				if (groupableProtein.getAccession().equals(acc)) {
					return groupableProtein.getEvidence();
				}
			}
		}
		return null;
	}

	private void writeHeaders(FileWriter fw, CellType celltype, Wash wash) throws IOException {

		// write the header
		fw.write("NUM\t" + "VALID\t" + "ACC\t" + "Gene\t" + "protein description\t" + "PROTEIN_EVIDENCE\t" + "SCORE\t"
				+ "Known NE\t" + "Transmembrane region\t" + "Total SPC in " + celltype + "/" + wash + "\t");

		final List<Experiment> experimentList = getAllExperiments();

		if (celltype != null) {
			for (final Experiment experiment : experimentList) {
				if (experiment.getCellType() == celltype) {
					if (wash == null || experiment.getWash() == wash) {
						// SPC replicates
						for (final CellCompartment cellCompartment : CellCompartment.values()) {
							final List<Replicate> replicates = experiment.getSortedReplicates();
							for (final Replicate replicate : replicates) {
								final Fractionation fractionation = replicate.getFractionation(cellCompartment);
								if (fractionation != null) {
									fw.write(fractionation.getName() + "_SPC\t");
								}
							}
						}
					}
				}
			}
			for (final Experiment experiment : experimentList) {
				if (experiment.getCellType() == celltype) {
					if (wash == null || experiment.getWash() == wash) {

						// NSAF replicates
						for (final CellCompartment cellCompartment : CellCompartment.values()) {
							final List<Replicate> replicates = experiment.getSortedReplicates();
							for (final Replicate replicate : replicates) {
								final Fractionation fractionation = replicate.getFractionation(cellCompartment);
								if (fractionation != null) {
									fw.write(fractionation.getName() + "_NSAF\t");
								}
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

			proteinAccs = new THashSet<String>();
			final List<Experiment> allExperiments = getAllExperiments();
			for (final Experiment experiment : allExperiments) {
				if (experiment.getCellType() == cellType) {
					proteinAccs.addAll(experiment.getProteinAccs());
				}
			}
		}
		return proteinAccs;
	}

	public List<String> getAccs(String rawAcc) {
		final List<String> accs = new ArrayList<>();

		// remove the [evidence] at the end
		if (rawAcc.contains("[")) {
			rawAcc = rawAcc.substring(0, rawAcc.indexOf("["));
		}
		if (rawAcc.contains(",")) {
			final String[] split = rawAcc.split(",");
			for (final String acc : split) {
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

	public String getProteinNameString(String rawAcc, CellType cellType, Wash wash) throws IOException {

		String proteinName = "";
		final List<Boolean> validArray = filterAccessionsByEvidence(rawAcc, null, cellType, wash);
		int index = 0;
		for (final String acc : getAccs(rawAcc)) {
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

	public String getGeneNameString(String rawAcc, CellType cellType, Wash wash) throws IOException {

		final List<String> geneNames = new ArrayList<String>();
		final List<Boolean> validArray = filterAccessionsByEvidence(rawAcc, null, cellType, wash);
		int index = 0;
		for (final String acc : getAccs(rawAcc)) {
			if (!validArray.get(index++)) {
				continue;
			}
			final Map<String, Protein> annotatedProtein = Constants.upr.getAnnotatedProtein(acc);
			if (annotatedProtein.containsKey(acc)) {
				final Protein protein2 = annotatedProtein.get(acc);
				if (protein2 != null) {
					boolean onePrimary = false;
					if (protein2.getGenes() != null) {
						for (final Gene gene : protein2.getGenes()) {
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
		final List<Experiment> list = new ArrayList<>();
		list.addAll(experimentsA);
		list.addAll(experimentsM);
		list.addAll(experimentsU);
		return list;
	}

	private ProteinGroup getGroup(String groupKey, CellType cellType, Wash wash) throws IOException {
		if (!groupsByRawAcc.containsKey(groupKey)) {
			final List<ProteinGroup> proteinGroupsTotal = getProteinGroups(cellType, wash);
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

	private String getFilteredProteinEvidences(String rawAcc, CellType cellType, Wash wash) throws IOException {
		if (rawAcc.contains("Q99M74")) {
			log.warn(rawAcc);
		}

		final List<Boolean> validArray = filterAccessionsByEvidence(rawAcc, null, cellType, wash);
		final List<ProteinEvidence> evidences = new ArrayList<ProteinEvidence>();
		int index = 0;
		final ProteinGroup group = getGroup(rawAcc, cellType, wash);
		// group can be null
		if (group == null) {
			return "";
		}
		for (final String acc : getAccs(rawAcc)) {
			if (!validArray.get(index++)) {
				continue;
			}
			final ProteinEvidence evidence = getEvidence(group, acc);
			if (!evidences.contains(evidence)) {
				evidences.add(evidence);
			}
		}
		final StringBuilder sb = new StringBuilder();

		for (final ProteinEvidence evidence : evidences) {
			if (!"".equals(sb.toString())) {
				sb.append(Constants.SEPARATOR);
			}
			if (evidence != null) {
				sb.append(evidence);
			} else {

				for (final GroupableProtein groupableProtein : group) {
					log.warn(groupableProtein.getAccession());
				}
				sb.append("");
			}
		}
		return sb.toString();
	}

	private String getProteinEvidences(String rawAcc) throws IOException {

		final List<String> proteinAccs = getAccs(rawAcc);

		final ProteinGroup group = getGroup(rawAcc, null, null);
		final StringBuilder sb = new StringBuilder();
		if (group != null) {
			for (final String proteinAcc : proteinAccs) {
				final ProteinEvidence evidence = getEvidence(group, proteinAcc);
				if (!"".equals(sb.toString())) {
					sb.append(Constants.SEPARATOR);
				}
				if (evidence != null) {
					sb.append(evidence);
				} else {
					log.warn(" I didnt find " + proteinAcc + "  in group  " + group);
					for (final GroupableProtein groupableProtein : group) {
						log.warn(groupableProtein.getAccession());
					}
					sb.append("");
				}
			}

		}
		return sb.toString();
	}

	private String getGeneName(List<String> geneNames) {
		final StringBuilder sb = new StringBuilder();

		final Set<String> set = new THashSet<String>();
		for (final String geneName : geneNames) {
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

	public List<Experiment> getExperiments(CellType cellType, Wash wash) {
		final List<Experiment> ret = new ArrayList<Experiment>();
		for (final Experiment experiment : getAllExperiments()) {
			if (cellType == null || experiment.getCellType() == cellType) {
				if (wash == null || wash == experiment.getWash()) {
					ret.add(experiment);
				}
			}
		}
		return ret;
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
			final ProteinGroup group = getGroup(rawAcc, null, null);
			if (group != null) {
				for (final GroupableProtein groupableProtein : group) {
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
