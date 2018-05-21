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
			Constants.writeCombinedDistribution = false;// UAM
			Constants.compareScores = false;
			UniprotProteinRetriever.enableCache = true;
			scoringFunction = new ScoringFunctionByNE_NSAF_Ratios(analyzer);
			// scoringFunction = new ScoringFunctionByNE_NSAF_Points(analyzer);
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

	private final List<Experiment> experimentsNoWash = new ArrayList<Experiment>();
	private final List<Experiment> experimentsUreaWash = new ArrayList<Experiment>();
	private final List<Experiment> experimentsCarbonateWash = new ArrayList<Experiment>();
	private final List<Experiment> experimentsPreWash = new ArrayList<Experiment>();

	public _4DNucleomeAnalyzerXI(String pass) throws IOException {
		super(pass);
		datasetsPathsFile = "Z:\\share\\Salva\\data\\4D_Nucleome\\Xi data\\Xi data May 2018\\IDcpr10_DTAselectTxt_FullPath_03122018Xi_cleaned.txt";
		outputFolder = new File(new File(datasetsPathsFile).getParent() + File.separator + "output");

	}

	private void loadDatasetsXi() throws IOException {
		log.info("Loading datasets");
		long t1 = System.currentTimeMillis();

		experimentsNoWash.clear();
		experimentsUreaWash.clear();
		experimentsCarbonateWash.clear();
		experimentsPreWash.clear();

		DataPaths dataPaths = new DataPaths(Constants.DATASET_PATHS_FILE);
		// U (N, Ne, C)

		// FDR 1%

		Pair<String, String> noWash = dataPaths.getXiFiles("mlNE_no wash_2D_rep1");
		Pair<String, String> cNE_Cw1 = dataPaths.getXiFiles("cNE_Cw1_2D_rep1");
		Pair<String, String> cCM_Cw1_LC1 = dataPaths.getXiFiles("cCM_Cw1_2D_rep1_LC1");
		Pair<String, String> cCM_Cw1_LC2 = dataPaths.getXiFiles("cCM_Cw1_2D_rep1_LC2");
		Pair<String, String> cNE_Uw1 = dataPaths.getXiFiles("cNE_Uw1_2D_rep1");
		Pair<String, String> cCM_Uw1 = dataPaths.getXiFiles("cCM_Uw1_2D_rep1");
		Pair<String, String> cNE_Uw2 = dataPaths.getXiFiles("cNE_Uw2_0216LT_rep2");
		Pair<String, String> cCM_Uw2 = dataPaths.getXiFiles("cCM_Uw2_0216LT_rep2");
		Pair<String, String> cCM_Cw1_3x = dataPaths.getXiFiles("cCM_Cw1_3x_dig2_0221LT");
		Pair<String, String> cCM_Uw1_3x = dataPaths.getXiFiles("cCM_Uw2_3x_dig2_0221LT");
		Pair<String, String> cNE_preCw1 = dataPaths.getXiFiles("20180316_01_cNE_preCw1_0303LT");
		Pair<String, String> cCM_preCw1 = dataPaths.getXiFiles("20180316_02_cCM_preCw1_0303LT");
		Pair<String, String> cNE_preUw1 = dataPaths.getXiFiles("20180316_03_cNE_preUw1_0303LT");
		Pair<String, String> cCM_preUw1 = dataPaths.getXiFiles("20180316_04_cCM_preUw1_0303LT");
		Pair<String, String> cNE_Uw3 = dataPaths.getXiFiles("20180406_01_cNE_Uw3_0403LT");
		Pair<String, String> cCM_Uw3 = dataPaths.getXiFiles("20180406_02_cCM_Uw3_0403LT");

		Experiment experimentnoWash = new Experiment("mlNE_no wash", Wash.NONE, CellType.LIVER);
		experimentnoWash.addReplicate(1, Wash.NONE, CellType.LIVER, CellCompartment.NE, getRemoteFile(noWash));
		experimentsNoWash.add(experimentnoWash);

		Experiment experimentCarbonateWash = new Experiment("CARBONATE_WASH", Wash.CARBONATE, CellType.C3H);
		experimentsCarbonateWash.add(experimentCarbonateWash);
		experimentCarbonateWash.addReplicate(1, Wash.CARBONATE, CellType.C3H, CellCompartment.NE,
				getRemoteFile(cNE_Cw1));
		experimentCarbonateWash.addReplicate(1, Wash.CARBONATE, CellType.C3H, CellCompartment.CM,
				getRemoteFile(cCM_Cw1_LC1));
		experimentCarbonateWash.addReplicate(2, Wash.CARBONATE, CellType.C3H, CellCompartment.CM,
				getRemoteFile(cCM_Cw1_LC2));
		Experiment experimentCarbonateWash3X = new Experiment("CARBONATE_WASH_3X", Wash.CARBONATE_3X, CellType.C3H);
		experimentCarbonateWash3X.addReplicate(1, Wash.CARBONATE_3X, CellType.C3H, CellCompartment.CM,
				getRemoteFile(cCM_Cw1_3x));
		experimentsCarbonateWash.add(experimentCarbonateWash3X);

		Experiment experimentUreaWash = new Experiment("UREA_WASH", Wash.UREA, CellType.C3H);
		experimentsUreaWash.add(experimentUreaWash);
		experimentUreaWash.addReplicate(1, Wash.UREA, CellType.C3H, CellCompartment.NE, getRemoteFile(cNE_Uw1));
		experimentUreaWash.addReplicate(2, Wash.UREA, CellType.C3H, CellCompartment.NE, getRemoteFile(cNE_Uw2));
		experimentUreaWash.addReplicate(3, Wash.UREA, CellType.C3H, CellCompartment.NE, getRemoteFile(cNE_Uw3));
		experimentUreaWash.addReplicate(1, Wash.UREA, CellType.C3H, CellCompartment.CM, getRemoteFile(cCM_Uw1));
		experimentUreaWash.addReplicate(2, Wash.UREA, CellType.C3H, CellCompartment.CM, getRemoteFile(cCM_Uw2));
		experimentUreaWash.addReplicate(3, Wash.UREA, CellType.C3H, CellCompartment.CM, getRemoteFile(cCM_Uw3));

		Experiment experimentUreaWash3x = new Experiment("UREA_WASH_3X", Wash.UREA_3X, CellType.C3H);
		experimentUreaWash3x.addReplicate(1, Wash.UREA_3X, CellType.C3H, CellCompartment.CM, getRemoteFile(cCM_Uw1_3x));
		experimentsUreaWash.add(experimentUreaWash3x);

		Experiment experimentPreCarbonateWash = new Experiment("PRE_CARBONATE_WASH", Wash.PREC, CellType.C3H);
		experimentsPreWash.add(experimentPreCarbonateWash);
		experimentPreCarbonateWash.addReplicate(1, Wash.PREC, CellType.C3H, CellCompartment.NE,
				getRemoteFile(cNE_preCw1));
		experimentPreCarbonateWash.addReplicate(1, Wash.PREC, CellType.C3H, CellCompartment.CM,
				getRemoteFile(cCM_preCw1));
		Experiment experimentPreUreaWash = new Experiment("PRE_UREA_WASH", Wash.PREU, CellType.C3H);
		experimentPreUreaWash.addReplicate(1, Wash.PREU, CellType.C3H, CellCompartment.NE, getRemoteFile(cNE_preUw1));
		experimentPreUreaWash.addReplicate(1, Wash.PREU, CellType.C3H, CellCompartment.CM, getRemoteFile(cCM_preUw1));

		long t2 = System.currentTimeMillis();
		log.info("It took " + DatesUtil.getDescriptiveTimeFromMillisecs(t2 - t1));

	}

	@Override
	public void run() throws IOException {
		// choose scoring function
		long t1 = System.currentTimeMillis();
		try {
			// load data
			loadDatasetsXi();
			// print scores for each cell type
			final CellType[] values = CellType.values();
			for (CellType cellType : values) {
				// restart field variables
				proteinAccs = null;
				proteinGroups = null;
				totalSPCs.clear();
				groupsByRawAcc.clear();
				groupableProteins.clear();
				for (Wash wash : Wash.values()) {
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
		List<Experiment> list = new ArrayList<Experiment>();
		list.addAll(experimentsNoWash);
		list.addAll(experimentsUreaWash);
		list.addAll(experimentsCarbonateWash);
		list.addAll(experimentsPreWash);
		return list;
	}

}
