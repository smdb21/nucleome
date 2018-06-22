package edu.scripps.yates.nucleome;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinRetriever;
import edu.scripps.yates.nucleome.model.CellCompartment;
import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.nucleome.model.Experiment;
import edu.scripps.yates.nucleome.model.Wash;
import edu.scripps.yates.utilities.dates.DatesUtil;
import edu.scripps.yates.utilities.util.Pair;

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
			Constants.DATASET_PATHS_FILE = datasetsPathsFile;
			Constants.MIN_TOTAL_SPC = 5;
			Constants.MAX_TEST_PROTEINS = 200000;
			Constants.writeCombinedDistribution = true;// UAM
			Constants.compareScores = false;
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

	private final List<Experiment> exps_uw1 = new ArrayList<Experiment>();
	private final List<Experiment> exps_uw2 = new ArrayList<Experiment>();
	private final List<Experiment> exps_uw3 = new ArrayList<Experiment>();
	private final List<Experiment> exps_uw4 = new ArrayList<Experiment>();

	private final List<Experiment> exps_cw1 = new ArrayList<Experiment>();
	private final List<Experiment> exps_cw3 = new ArrayList<Experiment>();
	private final List<Experiment> exps_preWC1 = new ArrayList<Experiment>();
	private final List<Experiment> exps_preWU1 = new ArrayList<Experiment>();
	private final List<Experiment> exps_preWU4 = new ArrayList<Experiment>();

	public _4DNucleomeAnalyzerXI(String pass) throws IOException {
		super(pass);
		datasetsPathsFile = "Z:\\share\\Salva\\data\\4D_Nucleome\\Xi data\\Xi data June 2018\\20180621_data_paths.txt";
		outputFolder = new File(new File(datasetsPathsFile).getParent() + File.separator + "output");

	}

	private void loadDatasetsXi() throws IOException {
		log.info("Loading datasets");
		final long t1 = System.currentTimeMillis();
		exps_uw1.clear();
		exps_uw2.clear();
		exps_uw3.clear();
		exps_cw1.clear();
		exps_cw3.clear();
		exps_preWC1.clear();
		exps_preWU1.clear();
		exps_preWU4.clear();
		final DataPaths dataPaths = new DataPaths(Constants.DATASET_PATHS_FILE);
		// U (N, Ne, C)

		// FDR 1%

		// final Pair<String, String> noWash = dataPaths.getXiFiles("mlNE_no
		// wash_2D_rep1");
		// Cw1
		final Pair<String, String> NE_C1 = dataPaths.getXiFiles("NE_C1");
		final Pair<String, String> CM_C1_3x = dataPaths.getXiFiles("CM_C1_3x");
		final Experiment cw1 = new Experiment("Cw1", Wash.CW1, CellType.C3H);
		exps_cw1.add(cw1);
		cw1.addReplicate(1, Wash.CW1, CellType.C3H, CellCompartment.NE, getRemoteFile(NE_C1));
		cw1.addReplicate(1, Wash.CW1, CellType.C3H, CellCompartment.CM, getRemoteFile(CM_C1_3x));
		// Cw3
		final Pair<String, String> CM_C3 = dataPaths.getXiFiles("CM_C3");
		final Pair<String, String> NE_C3 = dataPaths.getXiFiles("NE_C3");
		final Experiment cw3 = new Experiment("Cw3", Wash.CW3, CellType.C3H);
		exps_cw3.add(cw3);
		cw3.addReplicate(1, Wash.CW3, CellType.C3H, CellCompartment.CM, getRemoteFile(CM_C3));
		cw3.addReplicate(1, Wash.CW3, CellType.C3H, CellCompartment.NE, getRemoteFile(NE_C3));
		// preWC1
		final Pair<String, String> CM_preC1 = dataPaths.getXiFiles("CM_preC1");
		final Pair<String, String> NE_preC1 = dataPaths.getXiFiles("NE_preC1");
		final Experiment preWC1 = new Experiment("preWC1", Wash.PREWC1, CellType.C3H);
		preWC1.addReplicate(1, Wash.PREWC1, CellType.C3H, CellCompartment.CM, getRemoteFile(CM_preC1));
		preWC1.addReplicate(1, Wash.PREWC1, CellType.C3H, CellCompartment.NE, getRemoteFile(NE_preC1));
		exps_preWC1.add(preWC1);
		// preWU1
		final Pair<String, String> CM_preU1 = dataPaths.getXiFiles("CM_preU1");
		final Pair<String, String> NE_preU1 = dataPaths.getXiFiles("NE_preU1");
		final Experiment preWU1 = new Experiment("preWU1", Wash.preWU1, CellType.C3H);
		exps_preWU1.add(preWU1);
		preWU1.addReplicate(1, Wash.preWU1, CellType.C3H, CellCompartment.CM, getRemoteFile(CM_preU1));
		preWU1.addReplicate(1, Wash.preWU1, CellType.C3H, CellCompartment.NE, getRemoteFile(NE_preU1));
		// preWU4
		final Pair<String, String> CM_preU4 = dataPaths.getXiFiles("CM_preU4");
		final Pair<String, String> NE_preU4 = dataPaths.getXiFiles("NE_preU4");
		final Experiment preWU4 = new Experiment("preWU4", Wash.PREWU4, CellType.C3H);
		exps_preWU4.add(preWU4);
		preWU4.addReplicate(1, Wash.PREWU4, CellType.C3H, CellCompartment.CM, getRemoteFile(CM_preU4));
		preWU4.addReplicate(1, Wash.PREWU4, CellType.C3H, CellCompartment.NE, getRemoteFile(NE_preU4));
		// Uw1
		final Pair<String, String> CM_U1_3x_flush = dataPaths.getXiFiles("CM_U1_3x_flush");
		final Pair<String, String> NE_U1 = dataPaths.getXiFiles("NE_U1");
		final Experiment uw1 = new Experiment("Uw1", Wash.UW1, CellType.C3H);
		uw1.addReplicate(1, Wash.UW1, CellType.C3H, CellCompartment.CM, getRemoteFile(CM_U1_3x_flush));
		uw1.addReplicate(1, Wash.UW1, CellType.C3H, CellCompartment.NE, getRemoteFile(NE_U1));
		exps_uw1.add(uw1);
		// Uw2
		final Pair<String, String> CM_U2 = dataPaths.getXiFiles("CM_U2");
		final Pair<String, String> NE_U2 = dataPaths.getXiFiles("NE_U2");
		final Experiment uw2 = new Experiment("Uw2", Wash.UW2, CellType.C3H);
		uw2.addReplicate(1, Wash.UW2, CellType.C3H, CellCompartment.CM, getRemoteFile(CM_U2));
		uw2.addReplicate(1, Wash.UW2, CellType.C3H, CellCompartment.NE, getRemoteFile(NE_U2));
		exps_uw2.add(uw2);
		// Uw3
		final Pair<String, String> CM_U3 = dataPaths.getXiFiles("CM_U3");
		final Pair<String, String> NE_U3 = dataPaths.getXiFiles("NE_U3");
		final Experiment uw3 = new Experiment("Uw3", Wash.UW3, CellType.C3H);
		uw3.addReplicate(1, Wash.UW3, CellType.C3H, CellCompartment.CM, getRemoteFile(CM_U3));
		uw3.addReplicate(1, Wash.UW3, CellType.C3H, CellCompartment.NE, getRemoteFile(NE_U3));
		exps_uw3.add(uw3);
		// Uw4
		final Pair<String, String> CM_U4 = dataPaths.getXiFiles("CM_U4");
		final Pair<String, String> NE_U4 = dataPaths.getXiFiles("NE_U4");
		final Experiment uw4 = new Experiment("Uw4", Wash.UW4, CellType.C3H);
		uw4.addReplicate(1, Wash.UW4, CellType.C3H, CellCompartment.CM, getRemoteFile(CM_U4));
		uw4.addReplicate(1, Wash.UW4, CellType.C3H, CellCompartment.NE, getRemoteFile(NE_U4));
		exps_uw4.add(uw4);

		final long t2 = System.currentTimeMillis();
		log.info("It took " + DatesUtil.getDescriptiveTimeFromMillisecs(t2 - t1));

	}

	@Override
	public void run() throws IOException {
		// choose scoring function
		final long t1 = System.currentTimeMillis();
		try {
			// load data
			loadDatasetsXi();
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
					if (cellType != CellType.C3H) {
						continue;
					}

					// annotate proteins with uniprot
					annotateProteins(cellType);
					writeScoreDistributions(cellType, wash);

					// writeScoreDistributions(cellType, DataType.NSAF, true &&
					// Constants.printRatios);
					// writeScoreDistributions(cellType, DataType.PEPC, false &&
					// Constants.printRatios);
				}
			}
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
				writeScoreDistributions(null, null);
			}
			// writeScoreDistributions(null, DataType.PEPC, false &&
			// Constants.printRatios);

		} finally {
			log.info("It took " + DatesUtil.getDescriptiveTimeFromMillisecs(System.currentTimeMillis() - t1));
		}
	}

	@Override
	public List<Experiment> getAllExperiments() {
		final List<Experiment> list = new ArrayList<Experiment>();
		list.addAll(exps_uw1);
		list.addAll(exps_uw2);
		list.addAll(exps_uw3);
		list.addAll(exps_uw4);
		list.addAll(exps_cw1);
		list.addAll(exps_cw3);
		list.addAll(exps_preWC1);
		list.addAll(exps_preWU1);
		list.addAll(exps_preWU4);
		return list;
	}

}
