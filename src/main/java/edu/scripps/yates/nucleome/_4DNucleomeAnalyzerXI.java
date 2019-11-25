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

import org.apache.commons.lang.time.DateFormatUtils;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinRetriever;
import edu.scripps.yates.nucleome.model.CellCompartment;
import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.nucleome.model.ControlNE;
import edu.scripps.yates.nucleome.model.Experiment;
import edu.scripps.yates.nucleome.model.Fractionation;
import edu.scripps.yates.nucleome.model.Replicate;
import edu.scripps.yates.nucleome.model.Wash;
import edu.scripps.yates.nucleome.model.WashGroup;
import edu.scripps.yates.utilities.dates.DatesUtil;
import edu.scripps.yates.utilities.grouping.ProteinEvidence;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import smile.math.Histogram;

public class _4DNucleomeAnalyzerXI extends _4DNucleomeAnalyzer {
	private final static Logger log = Logger.getLogger(_4DNucleomeAnalyzerXI.class);

	public static void main(String[] args) {
		_4DNucleomeAnalyzerXI analyzer;
		try {

			analyzer = new _4DNucleomeAnalyzerXI(args[0]);

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
			Constants.MIN_TOTAL_SPC = 5;
			Constants.MAX_TEST_PROTEINS = 200000;
			Constants.writeCombinedDistribution = true;// UAM
			Constants.compareScores = false;
			UniprotProteinRetriever.enableCache = true;
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

	//
	private final List<Experiment> exps_u1 = new ArrayList<Experiment>();
	private final List<Experiment> exps_u2 = new ArrayList<Experiment>();
	private final List<Experiment> exps_u3 = new ArrayList<Experiment>();
	private final List<Experiment> exps_u4 = new ArrayList<Experiment>();

	//
	private final List<Experiment> exps_c1 = new ArrayList<Experiment>();
	private final List<Experiment> exps_c3 = new ArrayList<Experiment>();
	// private final List<Experiment> exps_c4 = new ArrayList<Experiment>();
	//
	private final List<Experiment> exps_preC1 = new ArrayList<Experiment>();
	//
	private final List<Experiment> exps_preU1 = new ArrayList<Experiment>();
	private final List<Experiment> exps_preU2 = new ArrayList<Experiment>();
	private final List<Experiment> exps_preU3 = new ArrayList<Experiment>();
	private final List<Experiment> exps_preU4 = new ArrayList<Experiment>();

	//
	// private final List<Experiment> exps_wc = new ArrayList<Experiment>();
	private DataPaths dataPaths;
	private final CellType cellType = CellType.C3H;
	private List<ProteinGroup> groups;
	private ArrayList<String> sortedAccsBySPC;
	private _4DNucleomeAnalyzer sabysDataAnalyzer;

	public _4DNucleomeAnalyzerXI(String pass) throws IOException {
		super(pass);
		setDatasetsPathsFile("Z:\\share\\Salva\\data\\4D_Nucleome\\Xi data\\Xi data Feb 2019\\20190201_data_paths.txt");

	}

	private Map<String, CalculatedScores> loadDatasetsXi() throws IOException {
		final Map<String, CalculatedScores> ret = new THashMap<String, CalculatedScores>();
		log.info("Loading datasets");
		final long t1 = System.currentTimeMillis();
		exps_u1.clear();
		exps_u2.clear();
		exps_u3.clear();
		exps_u4.clear();
		exps_c1.clear();
		exps_c3.clear();
		// exps_c4.clear();
		exps_preC1.clear();
		exps_preU1.clear();
		exps_preU3.clear();
		exps_preU4.clear();
		// exps_wc.clear();
		dataPaths = new DataPaths(getDatasetsPathsFile());
		// U (N, Ne, C)

		// FDR 1%

		// U1
		final Experiment u1 = new Experiment("u1", Wash.U1, cellType);
		exps_u1.add(u1);
		addReplicate(u1, "U1-CML");
		addReplicate(u1, "U1-NE");
		// U2
		final Experiment u2 = new Experiment("u2", Wash.U2, cellType);
		exps_u2.add(u2);
		addReplicate(u2, "U2-CML");
		addReplicate(u2, "U2-CMH new");
		addReplicate(u2, "U2-CMH old");
		addReplicate(u2, "U2-NE");
		// U3
		final Experiment u3 = new Experiment("u3", Wash.U3, cellType);
		exps_u3.add(u3);
		addReplicate(u3, "U3-CML");
		addReplicate(u3, "U3-CMH new");
		addReplicate(u3, "U3-CMH old");
		addReplicate(u3, "U3-NE");
		// U4
		final Experiment u4 = new Experiment("u4", Wash.U4, cellType);
		exps_u4.add(u4);
		addReplicate(u4, "U4-CML");
		addReplicate(u4, "U4-CMH new");
		addReplicate(u4, "U4-CMH old");
		addReplicate(u4, "U4-NE");

		// C1
		final Experiment c1 = new Experiment("c1", Wash.C1, cellType);
		exps_c1.add(c1);
		addReplicate(c1, "C1-CML");
		addReplicate(c1, "C1-NE");
		// C3
		final Experiment c3 = new Experiment("c3", Wash.C3, cellType);
		exps_c3.add(c3);
		addReplicate(c3, "C3-CML");
		addReplicate(c3, "C3-NE");
		// preC1
		final Experiment preC1 = new Experiment("preC1", Wash.preC1, cellType);
		exps_preC1.add(preC1);
		addReplicate(preC1, "preC1-CML");
		addReplicate(preC1, "preC1-NE");
		// preU1
		final Experiment preU1 = new Experiment("preU1", Wash.preU1, cellType);
		exps_preU1.add(preU1);
		addReplicate(preU1, "preU1-CML");
		addReplicate(preU1, "preU1-NE");
		// preU2
		final Experiment preU2 = new Experiment("preU2", Wash.preU2, cellType);
		exps_preU2.add(preU2);
		addReplicate(preU2, "preU2-CMH new");
		addReplicate(preU2, "preU2-CMH old");
		// preU3
		final Experiment preU3 = new Experiment("preU3", Wash.preU3, cellType);
		exps_preU3.add(preU3);
		addReplicate(preU3, "preU3-CMH new");
		addReplicate(preU3, "preU3-CMH old");
		// preU4
		final Experiment preU4 = new Experiment("preU4", Wash.preU4, cellType);
		exps_preU4.add(preU4);
		addReplicate(preU4, "preU4-CML");
		addReplicate(preU4, "preU4-CMH new");
		addReplicate(preU4, "preU4-CMH old");
		addReplicate(preU4, "preU4-NE");

		// WC
		// final Experiment wc = new Experiment("wc", Wash.WC, cellType);
		// exps_wc.add(wc);
		// addReplicate(wc, "WC-CML");
		// addReplicate(wc, "WC-NE");

		final long t2 = System.currentTimeMillis();
		log.info("It took " + DatesUtil.getDescriptiveTimeFromMillisecs(t2 - t1));

		// annotate proteins with uniprot
		annotateProteins(cellType);

		groups = getProteinGroups(cellType, null);
		for (final ProteinGroup proteinGroup : groups) {
			final String rawAcc = proteinGroup.getKey();
			final String acc = getAccessionStringByEvidence(rawAcc, proteinGroup, cellType, null);
			final CalculatedScores scores = new CalculatedScores(acc);
			ret.put(acc, scores);
			for (final String individualAcc : getAccs(rawAcc)) {
				if (ret.containsKey(individualAcc) && !individualAcc.equals(acc)) {
					log.info("asdf");
				}
				ret.put(individualAcc, scores);

			}
		}

		return ret;
	}

	private void addReplicate(Experiment experiment, String name) throws IOException {
		addReplicate(experiment, name, 1);
	}

	private void addReplicate(Experiment experiment, String name, int replicateNumber) throws IOException {
		final Pair<String, String> xiFile = dataPaths.getXiFiles(name);
		if (xiFile != null) {
			final File remoteFile = getRemoteFile(xiFile);
			Replicate replicate = experiment.getReplicateMap().get(replicateNumber);
			if (replicate == null) {
				replicate = new Replicate(experiment.getName(), replicateNumber, getWashFromName(name), cellType);
			}
			replicate.setFraction(getCellCompartmentFromName(name), getWashFromName(name),
					new RemoteSSHFileReference(remoteFile));
			experiment.addReplicate(replicate);
		} else {
			log.warn("File for " + name + " is not found");
		}
	}

	private Wash getWashFromName(String name) {
		return Wash.valueOf(name.split("-")[0]);
	}

	private CellCompartment getCellCompartmentFromName(String name) {
		final String tmp = name.split("-")[1];
		if (tmp.toLowerCase().contains("new")) {
			return CellCompartment.CMH_NEW;
		} else if (tmp.toLowerCase().contains("old")) {
			return CellCompartment.CMH_OLD;
		} else {
			return CellCompartment.valueOf(tmp);
		}
	}

	@Override
	public void run() throws IOException {
		// choose scoring function
		final long t1 = System.currentTimeMillis();
		try {
			// load data
			final Map<String, CalculatedScores> scores = loadDatasetsXi();

			for (final Wash wash : Wash.values()) {
				// restart field variables
				proteinAccs = null;
				totalSPCs.clear();
				groupsByRawAcc.clear();
				groupableProteins.clear();

				calculateIndividualScoresWithNoCMHFraction(cellType, wash, CellCompartment.NE, scores, groups);
				if (containsCMHFraction(wash)) {
					calculateIndividualScoresWithCMHFraction(cellType, wash, CellCompartment.NE, scores, groups);
				}
			}
			for (final WashGroup washGroup : WashGroup.values()) {
				calculateCompositeScore(washGroup, scores, groups);
			}

			calculateFoldChanges(scores, groups);
			calculateNuclearContentColumnsFromSabysDataset(scores, groups);
			writeScoreTableToFile(scores);
			writeScoresOnFileToCluster();
			// print scores for all celltypes together
			// writeScoreDistributions(null, DataType.NSAF, true &&
			// Constants.printRatios);
			// restart field variables
			proteinAccs = null;
			proteinGroups = null;
			totalSPCs.clear();
			groupsByRawAcc.clear();
			groupableProteins.clear();
			if (Constants.writeCombinedDistribution) {
				// writeScoreDistributions(null, null, CellCompartment.NE);
			}
			// writeScoreDistributions(null, DataType.PEPC, false &&
			// Constants.printRatios);

		} finally {
			log.info("It took " + DatesUtil.getDescriptiveTimeFromMillisecs(System.currentTimeMillis() - t1));
		}
	}

	private void writeScoresOnFileToCluster() throws IOException {
		final String formatedDate = DateFormatUtils.format(new Date(), "yyyy-MM-dd_HH-mm");
		final String pathname = getOutputFolder().getAbsolutePath() + File.separator + formatedDate + "_" + cellType
				+ "_ScoresToClusterTable_" + scoringFunction.getName() + ".txt";

		final File scoreFileOutput = new File(pathname);
		if (!scoreFileOutput.getParentFile().exists()) {
			log.info("Creating file " + scoreFileOutput.getAbsolutePath());
			scoreFileOutput.getParentFile().mkdirs();
		}
		if (scoreFileOutput.exists()) {
			log.info("Overriding file " + scoreFileOutput.getAbsolutePath());
		}
		final Set<Wash> washes = new THashSet<Wash>();
		washes.add(Wash.U2);
		washes.add(Wash.U3);
		washes.add(Wash.U4);
		final Set<CellCompartment> fractionTypes = new THashSet<CellCompartment>();
		fractionTypes.add(CellCompartment.CML);
		fractionTypes.add(CellCompartment.CMH_NEW);
		fractionTypes.add(CellCompartment.CMH_OLD);
		fractionTypes.add(CellCompartment.NE);

		FileWriter fw = null;
		try {
			final ProgressCounter counter = new ProgressCounter(sortedAccsBySPC.size(),
					ProgressPrintingType.PERCENTAGE_STEPS, 0, true);
			fw = new FileWriter(scoreFileOutput);
			fw.write("Num\tACC\tGene\tProtein Name\tTotal SPC\t");
			fw.write("CM\tCMH new\tCMH old\tNE\n");
			int num = 1;
			for (final String proteinGroupKey : sortedAccsBySPC) {
				counter.increment();
				final String printIfNecessary = counter.printIfNecessary();
				if (!"".equals(printIfNecessary)) {
					log.info(printIfNecessary);
				}
				final String rawAccString = proteinGroupKey;
				String geneNameString = getGeneNameString(rawAccString, cellType, null);
				if ("".equals(geneNameString)) {
					geneNameString = rawAccString;
				}
				final String proteinNameString = getProteinNameString(rawAccString, cellType, null);
				final int totalSPC = getTotalSPC(rawAccString, cellType, null);
				fw.write(num++ + "\t" + proteinGroupKey + "\t" + geneNameString + "\t" + proteinNameString + "\t"
						+ totalSPC + "\t");
				// CM score
				final List<String> proteinAccessions = getAccs(rawAccString);
				final double scoreCML = scoringFunction.getScore(proteinAccessions, cellType, washes,
						CellCompartment.CML, fractionTypes);
				fw.write(parseInfinities(scoreCML) + "\t");
				final double scoreCMHnew = scoringFunction.getScore(proteinAccessions, cellType, washes,
						CellCompartment.CMH_NEW, fractionTypes);
				fw.write(parseInfinities(scoreCMHnew) + "\t");
				final double scoreCMHOld = scoringFunction.getScore(proteinAccessions, cellType, washes,
						CellCompartment.CMH_OLD, fractionTypes);
				fw.write(parseInfinities(scoreCMHOld) + "\t");
				final double scoreNE = scoringFunction.getScore(proteinAccessions, cellType, washes, CellCompartment.NE,
						fractionTypes);
				fw.write(parseInfinities(scoreNE) + "\n");
			}
		} finally {
			fw.close();
		}
	}

	private void writeScoreTableToFile(Map<String, CalculatedScores> scores) throws IOException {
		final String formatedDate = DateFormatUtils.format(new Date(), "yyyy-MM-dd_HH-mm");
		final String pathname = getOutputFolder().getAbsolutePath() + File.separator + formatedDate + "_" + cellType
				+ "_TotalTable_" + scoringFunction.getName() + ".txt";

		final File scoreFileOutput = new File(pathname);
		if (!scoreFileOutput.getParentFile().exists()) {
			log.info("Creating file " + scoreFileOutput.getAbsolutePath());
			scoreFileOutput.getParentFile().mkdirs();
		}
		if (scoreFileOutput.exists()) {
			log.info("Overriding file " + scoreFileOutput.getAbsolutePath());
		}
		sortedAccsBySPC = new ArrayList<String>();
		sortedAccsBySPC.addAll(scores.keySet());
		log.info("Sorting by SPC");
		Collections.sort(sortedAccsBySPC, new Comparator<String>() {

			@Override
			public int compare(String o1, String o2) {
				try {
					final CalculatedScores scores1 = scores.get(o1);
					final CalculatedScores scores2 = scores.get(o2);
					final int totalSPC1 = getTotalSPC(scores1.getProteinKey(), cellType, null);
					final int totalSPC2 = getTotalSPC(scores2.getProteinKey(), cellType, null);
					return Integer.compare(totalSPC2, totalSPC1);
				} catch (final Exception e) {
					e.printStackTrace();
				}
				return 0;
			}
		});

		FileWriter fw = null;
		try {
			final ProgressCounter counter = new ProgressCounter(scores.size(), ProgressPrintingType.PERCENTAGE_STEPS, 0,
					true);
			fw = new FileWriter(scoreFileOutput);
			writeHeaders(fw, cellType);
			int num = 1;
			for (final String proteinGroupKey : sortedAccsBySPC) {
				counter.increment();
				final String printIfNecessary = counter.printIfNecessary();
				if (!"".equals(printIfNecessary)) {
					log.info(printIfNecessary);
				}
				final String rawAccString = proteinGroupKey;
				String geneNameString = getGeneNameString(rawAccString, cellType, null);
				if ("".equals(geneNameString)) {
					geneNameString = rawAccString;
				}
				final String proteinNameString = getProteinNameString(rawAccString, cellType, null);
				final CalculatedScores calculatedScores = scores.get(proteinGroupKey);
				final int totalSPC = getTotalSPC(rawAccString, cellType, null);
				final String valid = isValid(rawAccString, totalSPC) ? "VALID" : "FILTERED";
				final String proteinEvidence = getFilteredProteinEvidences(rawAccString, cellType, null);
				// write the header
				fw.write(num++ + "\t" + valid + "\t" + proteinGroupKey + "\t" + geneNameString + "\t"
						+ proteinNameString + "\t" + proteinEvidence + "\t");
				// individual scores
				for (final Wash wash : Wash.values()) {
					fw.write(parseInfinities(calculatedScores.getIndividualScoreWithNoCMHFraction(wash)) + "\t");
				}
				// scores with both CMH

				fw.write(parseInfinities(calculatedScores.getIndividualScoreWithCMHFraction(Wash.U2)) + "\t"
						+ parseInfinities(calculatedScores.getIndividualScoreWithCMHFraction(Wash.U3)) + "\t"
						+ parseInfinities(calculatedScores.getIndividualScoreWithCMHFraction(Wash.U4)) + "\t"
						+ parseInfinities(calculatedScores.getIndividualScoreWithCMHFraction(Wash.preU4)) + "\t");

				// composite scores
				for (final WashGroup washGroup : WashGroup.values()) {
					fw.write(calculatedScores.getGroupedScore(washGroup) + "\t");
				}

				// fold changes

				fw.write(parseInfinities(calculatedScores.getNEFoldChange(Wash.U1, Wash.preU1)) + "\t"
						+ parseInfinities(calculatedScores.getNEFoldChangePValue(Wash.U1, Wash.preU1)) + "\t"
						+ parseInfinities(calculatedScores.getNEFoldChange(Wash.C1, Wash.preC1)) + "\t"
						+ parseInfinities(calculatedScores.getNEFoldChangePValue(Wash.C1, Wash.preC1)) + "\t"
						+ parseInfinities(calculatedScores.getNEFoldChange(Wash.U4, Wash.preU4)) + "\t"
						+ parseInfinities(calculatedScores.getNEFoldChangePValue(Wash.U4, Wash.preU4)) + "\t"
						+ parseInfinities(calculatedScores.getNEFoldChange(Wash.C3, Wash.preU4)) + "\t"
						+ parseInfinities(calculatedScores.getNEFoldChangePValue(Wash.C3, Wash.preU4)) + "\t");
				//
				// new information requested from Li-Chun about Saby's dataset
				fw.write(calculatedScores.getTotalNuclearContentFractionSPCInUInSabyDataset() + "\t");
				fw.write(parseInfinities(calculatedScores.getNuclearContentEnrichmentScoreInUInSabyDataset()) + "\t");
				//
				for (final Wash wash : Wash.values()) {
					final int totalSPC2 = getTotalSPC(rawAccString, cellType, wash);
					fw.write(totalSPC2 + "\t");
				}
				//
				fw.write(ControlNE.isControl(rawAccString) + "\t" + getTransmembraneRegion(rawAccString, cellType, null)
						+ "\t" + totalSPC + "\t");

				final List<Experiment> experimentList = getAllExperiments();
				final List<String> filteredAccessions = getAccs(rawAccString);

				// SPC replicates
				for (final Experiment experiment : experimentList) {
					for (final CellCompartment cellCompartment : CellCompartment.values()) {
						final List<Replicate> replicates = experiment.getSortedReplicates();
						for (final Replicate replicate : replicates) {
							final Fractionation fractionation = replicate.getFractionation(cellCompartment);
							if (fractionation != null) {
								final int repSPC = fractionation.getSpectralCount(filteredAccessions, true);
								fw.write(repSPC + "\t");
							}
						}
					}
				}
				// NSAF replicates
				for (final Experiment experiment : experimentList) {
					for (final CellCompartment cellCompartment : CellCompartment.values()) {
						final List<Replicate> replicates = experiment.getSortedReplicates();
						for (final Replicate replicate : replicates) {
							final Fractionation fractionation = replicate.getFractionation(cellCompartment);
							if (fractionation != null) {
								final double repNSAF = fractionation.getAverageNSAF(filteredAccessions, true);
								fw.write(repNSAF + "\t");
							}
						}
					}
				}

				fw.write("\n");
			}
		} finally {
			fw.close();
		}

	}

	private String parseInfinities(double neFoldChange) {
		if (Double.POSITIVE_INFINITY == neFoldChange) {
			return ">1000";
		}
		if (Double.NEGATIVE_INFINITY == neFoldChange) {
			return "<1000";
		}
		if (Double.isNaN(neFoldChange)) {
			return "-";
		}
		return String.valueOf(neFoldChange);
	}

	private void writeHeaders(FileWriter fw, CellType celltype) throws IOException {

		// write the header
		fw.write("NUM\t" + "VALID\t" + "ACC\t" + "Gene\t" + "protein description\t" + "PROTEIN_EVIDENCE\t");
		// individual scores
		for (final Wash wash : Wash.values()) {
			fw.write("NE_Score_" + wash + "\t");
		}
		// scores with both CMH
		fw.write("NE_Score_U2_with_CMH\tNE_Score_U3_with_CMH\tNE_Score_U4_with_CMH\tNE_Score_PreU4_with_CMH\t");

		// composite scores
		for (final WashGroup washGroup : WashGroup.values()) {
			fw.write("NeScore_" + washGroup + "\t");
		}

		// fold changes
		fw.write("fold change U1/PreU1\t" + "p-value U1/PreU1\t" + "fold change C1/PreC1\t" + "p-value C1/PreC1\t"
				+ "fold change U4/PreU4\t" + "p-value U4/PreU4\t" + "fold change C3/PreU4\t" + "p-value C3/PreU4\t");
		// new data requested by Li-Chun
		fw.write("Total CN_SPC in U (Saby)\t");
		fw.write("NC score in U (Saby)\t");
		//
		for (final Wash wash : Wash.values()) {
			fw.write(wash + "_SPC\t");
		}
		//

		//
		fw.write("Known NE\t" + "Transmembrane region\tTotal SPC\t");

		final List<Experiment> experimentList = getAllExperiments();
		// SPC replicates
		for (final Experiment experiment : experimentList) {
			if (experiment.getCellType() == celltype) {
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
		// NSAF replicates
		for (final Experiment experiment : experimentList) {
			if (experiment.getCellType() == celltype) {
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

		fw.write("\n");

	}

	private void calculateNuclearContentColumnsFromSabysDataset(Map<String, CalculatedScores> scores,
			List<ProteinGroup> groups2) throws IOException {
		final _4DNucleomeAnalyzer sabysData = getSabysDataAnalyzer();
		sabysData.loadDatasets();
		for (final ProteinGroup proteinGroup : groups2) {
			if (proteinGroup.getEvidence() == ProteinEvidence.NONCONCLUSIVE) {
				continue;
			}
			final String rawAcc = proteinGroup.getKey();

			final List<String> filteredAcessions = new ArrayList<String>();
			final String filteredAcc = getAccessionStringByEvidence(rawAcc, proteinGroup, cellType, null);
			if (filteredAcc.contains(",")) {
				final String[] split = filteredAcc.split(",");
				for (final String acc : split) {
					filteredAcessions.add(acc);
				}
			} else {
				filteredAcessions.add(filteredAcc);
			}

			final List<Experiment> uExperiments = sabysData.getExperimentsU();
			int totalSPC = 0;

			for (final Experiment uExperiment : uExperiments) {
				final int spc = uExperiment.getSumSPC(filteredAcessions, CellCompartment.N, true);
				totalSPC += spc;
			}
			final double nuclearContentScore = sabysData.scoringFunction.getScore(filteredAcessions, CellType.U,
					CellCompartment.N);
			for (final String acc : filteredAcessions) {
				scores.get(acc).setTotalNuclearContentFractionSPCInUInSabyDataset(totalSPC);
				scores.get(acc).setNuclearContentEnrichmentScoreInUInSabyDataset(nuclearContentScore);
			}

		}
	}

	private _4DNucleomeAnalyzer getSabysDataAnalyzer() throws IOException {
		if (sabysDataAnalyzer == null) {
			sabysDataAnalyzer = new _4DNucleomeAnalyzer(pass);
		}
		return sabysDataAnalyzer;
	}

	private void calculateFoldChanges(Map<String, CalculatedScores> scores, List<ProteinGroup> groups2)
			throws IOException {
		final List<Pair<Wash, Wash>> listOfPairs = new ArrayList<Pair<Wash, Wash>>();

		listOfPairs.add(new Pair<Wash, Wash>(Wash.preU1, Wash.U1));
		listOfPairs.add(new Pair<Wash, Wash>(Wash.preC1, Wash.C1));
		listOfPairs.add(new Pair<Wash, Wash>(Wash.preU4, Wash.U4));
		listOfPairs.add(new Pair<Wash, Wash>(Wash.preU4, Wash.C3));
		for (final Pair<Wash, Wash> pair : listOfPairs) {
			final Wash washPreWash = pair.getFirstelement();
			final Wash washPostWash = pair.getSecondElement();
			final TDoubleArrayList ratioCollection = new TDoubleArrayList();
			for (final ProteinGroup proteinGroup : groups2) {
				if (proteinGroup.getEvidence() == ProteinEvidence.NONCONCLUSIVE) {
					continue;
				}
				final String rawAcc = proteinGroup.getKey();

				final List<String> filteredAcessions = new ArrayList<String>();
				final String filteredAcc = getAccessionStringByEvidence(rawAcc, proteinGroup, cellType, null);
				if (filteredAcc.contains(",")) {
					final String[] split = filteredAcc.split(",");
					for (final String acc : split) {
						filteredAcessions.add(acc);
					}
				} else {
					filteredAcessions.add(filteredAcc);
				}

				final double nsafPre = scoringFunction.getSumNSAFs(filteredAcessions, cellType, washPreWash,
						CellCompartment.NE);
				final double nsafPost = scoringFunction.getSumNSAFs(filteredAcessions, cellType, washPostWash,
						CellCompartment.NE);
				final double foldChange = nsafPost / nsafPre;

				final double log2ratio = Maths.log(foldChange, 2);
				if (Double.isFinite(log2ratio)) {
					ratioCollection.add(log2ratio);
				}
				scores.get(filteredAcc).addNEFoldChange(foldChange, washPostWash, washPreWash);

			}
			double[][] histogram = Histogram.histogram(ratioCollection.toArray());
			GaussianCurveFitter gaussianfitter = GaussianCurveFitter.create();
			WeightedObservedPoints obs = new WeightedObservedPoints();
			for (int j = 0; j < histogram[2].length; j++) {
				obs.add((histogram[0][j] + histogram[1][j]) / 2, histogram[2][j]);
			}
			double[] bestFit = gaussianfitter.fit(obs.toList());
			// final double gaussianNorm = bestFit[0];
			final double gaussianMean = bestFit[1];
			final double gaussianSigma = bestFit[2];
			// correct to center distribution on 0
			int index = 0;
			for (final double ratio : ratioCollection.toArray()) {
				ratioCollection.set(index++, ratio - gaussianMean);
			}
			histogram = Histogram.histogram(ratioCollection.toArray());
			gaussianfitter = GaussianCurveFitter.create();
			obs = new WeightedObservedPoints();
			for (int j = 0; j < histogram[2].length; j++) {
				obs.add((histogram[0][j] + histogram[1][j]) / 2, histogram[2][j]);
			}
			bestFit = gaussianfitter.fit(obs.toList());
			// final double gaussianNorm = bestFit[0];
			final double newGaussianMean = bestFit[1];
			final double newGaussianSigma = bestFit[2];
			final NormalDistribution normalDistribution = new NormalDistribution(gaussianMean, gaussianSigma);
			for (final ProteinGroup proteinGroup : groups2) {
				if (proteinGroup.getEvidence() == ProteinEvidence.NONCONCLUSIVE) {
					continue;
				}
				final String rawAcc = proteinGroup.getKey();

				final String filteredAcc = getAccessionStringByEvidence(rawAcc, proteinGroup, cellType, null);
				final double foldChange = scores.get(filteredAcc).getNEFoldChange(washPostWash, washPreWash);
				final double log2ratio = Maths.log(foldChange, 2);

				double pvalue = Double.NaN;
				if (Double.isInfinite(log2ratio)) {
					pvalue = 0;
				} else if (Double.isNaN(log2ratio)) {
					pvalue = Double.NaN;
				} else {
					pvalue = normalDistribution.cumulativeProbability(log2ratio - gaussianMean);
					// this p is the probability that a ratio is <= than ratio
					if (log2ratio > gaussianMean) {
						pvalue = 1 - pvalue; // this would be the probability
												// that a ratio is > than the
												// ratio
					}
				}
				scores.get(filteredAcc).addNEFoldChangePValue(pvalue, washPostWash, washPreWash);
			}
		}

	}

	private void calculateCompositeScore(WashGroup washGroup, Map<String, CalculatedScores> scores,
			List<ProteinGroup> groups2) throws IOException {

		final List<Wash> washes = Wash.getWashesWithAllFractionsFromWashGroup(washGroup);
		final Set<CellCompartment> fractionTypes = new THashSet<CellCompartment>();
		for (final CellCompartment cellCompartment : CellCompartment.getCellCompartmentsForXiDatasets()) {
			fractionTypes.add(cellCompartment);
		}

		final List<Pair<String, Double>> scores2 = calculateScoresFromGroups(cellType, washes, CellCompartment.NE,
				fractionTypes, groups2);
		for (final Pair<String, Double> pair : scores2) {
			final String groupKey = pair.getFirstelement();
			final double score = pair.getSecondElement();
			scores.get(groupKey).addGroupedScore(score, washGroup);
		}

	}

	private boolean containsCMHFraction(Wash wash) {
		switch (wash) {
		case U2:
		case U3:
		case U4:
		case preU4:

			return true;
		default:
			return false;
		}
	}

	protected void calculateIndividualScoresWithNoCMHFraction(CellType celltype, Wash wash,
			CellCompartment cellCompartmentToStudy, Map<String, CalculatedScores> scores,
			List<ProteinGroup> proteinGroups) throws IOException {

		final Set<CellCompartment> fractionTypes = new THashSet<CellCompartment>();
		fractionTypes.add(CellCompartment.CML);
		fractionTypes.add(CellCompartment.N);
		fractionTypes.add(CellCompartment.NE);
		final List<Pair<String, Double>> scores2 = calculateScoresFromGroups(celltype, wash, cellCompartmentToStudy,
				fractionTypes, proteinGroups);
		for (final Pair<String, Double> pair : scores2) {
			final String groupKey = pair.getFirstelement();
			final double score = pair.getSecondElement();
			if (scores.containsKey(groupKey)) {
				scores.get(groupKey).addIndividualScoreWithNoCMHFraction(score, wash);
			} else {
				throw new IllegalArgumentException(groupKey);
			}
		}

	}

	protected void calculateIndividualScoresWithCMHFraction(CellType celltype, Wash wash,
			CellCompartment cellCompartmentToStudy, Map<String, CalculatedScores> scores,
			List<ProteinGroup> proteinGroups) throws IOException {

		final Set<CellCompartment> fractionTypes = new THashSet<CellCompartment>();
		fractionTypes.add(CellCompartment.CML);
		fractionTypes.add(CellCompartment.N);
		fractionTypes.add(CellCompartment.NE);
		fractionTypes.add(CellCompartment.CMH_NEW);
		fractionTypes.add(CellCompartment.CMH_OLD);
		final List<Pair<String, Double>> scores2 = calculateScoresFromGroups(celltype, wash, cellCompartmentToStudy,
				fractionTypes, proteinGroups);
		for (final Pair<String, Double> pair : scores2) {
			final String groupKey = pair.getFirstelement();
			final double score = pair.getSecondElement();
			scores.get(groupKey).addIndividualScoreWithCMHFraction(score, wash);
		}

	}

	@Override
	public List<Experiment> getAllExperiments() {
		final List<Experiment> list = new ArrayList<Experiment>();
		list.addAll(exps_u1);
		list.addAll(exps_u2);
		list.addAll(exps_u3);
		list.addAll(exps_u4);
		list.addAll(exps_c1);
		list.addAll(exps_c3);

		// list.addAll(exps_c4);
		list.addAll(exps_preC1);
		list.addAll(exps_preU1);
		list.addAll(exps_preU4);
		list.addAll(exps_preU2);
		list.addAll(exps_preU3);

		// list.addAll(exps_wc);
		return list;
	}

}
