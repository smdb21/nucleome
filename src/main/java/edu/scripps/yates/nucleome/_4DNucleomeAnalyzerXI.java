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

public class _4DNucleomeAnalyzerXI extends _4DNucleomeAnalyzer {
	private final static Logger log = Logger.getLogger(_4DNucleomeAnalyzerXI.class);

	public static void main(String[] args) {
		_4DNucleomeAnalyzerXI analyzer;
		try {

			analyzer = new _4DNucleomeAnalyzerXI();

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

	public _4DNucleomeAnalyzerXI() throws IOException {
		datasetsPathsFile = "Z:\\share\\Salva\\data\\4D_Nucleome\\Xi data March18\\IDcpr10_DTAselectTxt_FullPath_03122018Xi.txt";
		outputFolder = new File(new File(datasetsPathsFile).getParent() + File.separator + "output");

	}

	private void loadDatasetsXi() throws IOException {
		log.info("Loading datasets");
		long t1 = System.currentTimeMillis();

		experimentsNoWash.clear();
		experimentsUreaWash.clear();
		experimentsCarbonateWash.clear();

		DataPaths dataPaths = new DataPaths(Constants.DATASET_PATHS_FILE);
		// U (N, Ne, C)

		// FDR 1%

		String noWash = dataPaths.getXiFiles("mlNE_no wash_2D_rep1");
		String cNE_Cw1 = dataPaths.getXiFiles("cNE_Cw1_2D_rep1");
		String cCM_Cw1_LC1 = dataPaths.getXiFiles("cCM_Cw1_2D_rep1_LC1");
		String cCM_Cw1_LC2 = dataPaths.getXiFiles("cCM_Cw1_2D_rep1_LC2");
		String cNE_Uw1 = dataPaths.getXiFiles("cNE_Uw1_2D_rep1");
		String cCM_Uw1 = dataPaths.getXiFiles("cCM_Uw1_2D_rep1");
		String cNE_Uw2 = dataPaths.getXiFiles("cNE_Uw2_0216LT_rep2");
		String cCM_Uw2 = dataPaths.getXiFiles("cCM_Uw2_0216LT_rep2");
		String cCM_Cw1_3x = dataPaths.getXiFiles("cCM_Cw1_3x_dig2_0221LT");
		String cCM_Uw1_3x = dataPaths.getXiFiles("cCM_Uw2_3x_dig2_0221LT");

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
		Experiment experimentCarbonateWash3X = new Experiment("CARBONATE_WASH_3X", Wash.CARBONATE, CellType.C3H);
		experimentCarbonateWash3X.addReplicate(1, Wash.CARBONATE, CellType.C3H, CellCompartment.CM,
				getRemoteFile(cCM_Cw1_3x));
		experimentsCarbonateWash.add(experimentCarbonateWash3X);

		Experiment experimentUreaWash = new Experiment("UREA_WASH", Wash.UREA, CellType.C3H);
		experimentsUreaWash.add(experimentUreaWash);
		experimentUreaWash.addReplicate(1, Wash.UREA, CellType.C3H, CellCompartment.NE, getRemoteFile(cNE_Uw1));
		experimentUreaWash.addReplicate(2, Wash.UREA, CellType.C3H, CellCompartment.NE, getRemoteFile(cNE_Uw2));
		experimentUreaWash.addReplicate(1, Wash.UREA, CellType.C3H, CellCompartment.CM, getRemoteFile(cCM_Uw1));
		experimentUreaWash.addReplicate(2, Wash.UREA, CellType.C3H, CellCompartment.CM, getRemoteFile(cCM_Uw2));

		Experiment experimentUreaWash3x = new Experiment("UREA_WASH_3X", Wash.UREA, CellType.C3H);
		experimentUreaWash3x.addReplicate(1, Wash.UREA, CellType.C3H, CellCompartment.CM, getRemoteFile(cCM_Uw1_3x));
		experimentsUreaWash.add(experimentUreaWash3x);

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

		return list;
	}

}
