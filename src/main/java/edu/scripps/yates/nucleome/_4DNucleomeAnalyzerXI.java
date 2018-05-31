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
		} catch (final IOException e) {
			e.printStackTrace();
			System.err.println("ERROR: " + e.getMessage());
		}
		System.exit(-1);
	}

	private final List<Experiment> exps_uw1 = new ArrayList<Experiment>();
	private final List<Experiment> exps_uw2 = new ArrayList<Experiment>();
	private final List<Experiment> exps_uw3 = new ArrayList<Experiment>();
	private final List<Experiment> exps_cw1 = new ArrayList<Experiment>();
	private final List<Experiment> exps_prew1 = new ArrayList<Experiment>();
	private final List<Experiment> exps_prew2 = new ArrayList<Experiment>();

	public _4DNucleomeAnalyzerXI(String pass) throws IOException {
		super(pass);
		datasetsPathsFile = "Z:\\share\\Salva\\data\\4D_Nucleome\\Xi data\\Xi data May 2018\\IDcpr10_DTAselectTxt_FullPath_03122018Xi_cleaned.txt";
		outputFolder = new File(new File(datasetsPathsFile).getParent() + File.separator + "output");

	}

	private void loadDatasetsXi() throws IOException {
		log.info("Loading datasets");
		final long t1 = System.currentTimeMillis();
		exps_uw1.clear();
		exps_uw2.clear();
		exps_uw3.clear();
		exps_cw1.clear();
		exps_prew1.clear();
		exps_prew2.clear();
		final DataPaths dataPaths = new DataPaths(Constants.DATASET_PATHS_FILE);
		// U (N, Ne, C)

		// FDR 1%

		// final Pair<String, String> noWash = dataPaths.getXiFiles("mlNE_no
		// wash_2D_rep1");
		final Pair<String, String> s3_Cw1_CM_3x = dataPaths.getXiFiles("#3_Cw1_CM-3x");
		final Pair<String, String> s3_Cw1_NE = dataPaths.getXiFiles("#3_Cw1_NE");
		final Pair<String, String> s3_preW1_CM = dataPaths.getXiFiles("#3_preW1_CM");
		final Pair<String, String> s3_preW1_NE = dataPaths.getXiFiles("#3_preW1_NE");
		final Pair<String, String> s4_preW2_CM = dataPaths.getXiFiles("#4_preW2_CM");
		final Pair<String, String> s4_preW2_NE = dataPaths.getXiFiles("#4_preW2_NE");
		final Pair<String, String> s4_Uw1_CM_3x = dataPaths.getXiFiles("#4_Uw1_CM-3x");
		final Pair<String, String> s4_Uw1_NE = dataPaths.getXiFiles("#4_Uw1_NE");
		final Pair<String, String> s6_Uw2_CM = dataPaths.getXiFiles("#6_Uw2_CM");
		final Pair<String, String> s6_Uw2_NE = dataPaths.getXiFiles("#6_Uw2_NE");

		final Pair<String, String> s12_Uw3_CM = dataPaths.getXiFiles("#12_Uw3_CM");
		final Pair<String, String> s12_Uw3_NE = dataPaths.getXiFiles("#12_Uw3_NE");

		final Experiment cw1 = new Experiment("Cw1", Wash.CW1, CellType.C3H);
		exps_cw1.add(cw1);
		cw1.addReplicate(1, Wash.CW1, CellType.C3H, CellCompartment.CM, getRemoteFile(s3_Cw1_CM_3x));
		cw1.addReplicate(1, Wash.CW1, CellType.C3H, CellCompartment.NE, getRemoteFile(s3_Cw1_NE));

		final Experiment preW1 = new Experiment("preW1", Wash.PREW1, CellType.C3H);
		preW1.addReplicate(1, Wash.PREW1, CellType.C3H, CellCompartment.CM, getRemoteFile(s3_preW1_CM));
		preW1.addReplicate(1, Wash.PREW1, CellType.C3H, CellCompartment.NE, getRemoteFile(s3_preW1_NE));
		exps_prew1.add(preW1);

		final Experiment preW2 = new Experiment("preW2", Wash.PREW2, CellType.C3H);
		exps_prew2.add(preW2);
		preW2.addReplicate(1, Wash.PREW2, CellType.C3H, CellCompartment.CM, getRemoteFile(s4_preW2_CM));
		preW2.addReplicate(1, Wash.PREW2, CellType.C3H, CellCompartment.NE, getRemoteFile(s4_preW2_NE));

		final Experiment uw1 = new Experiment("Uw1", Wash.UW1, CellType.C3H);
		uw1.addReplicate(1, Wash.UW1, CellType.C3H, CellCompartment.CM, getRemoteFile(s4_Uw1_CM_3x));
		uw1.addReplicate(1, Wash.UW1, CellType.C3H, CellCompartment.NE, getRemoteFile(s4_Uw1_NE));
		exps_uw1.add(uw1);

		final Experiment uw2 = new Experiment("Uw2", Wash.UW2, CellType.C3H);
		uw2.addReplicate(1, Wash.UW2, CellType.C3H, CellCompartment.CM, getRemoteFile(s6_Uw2_CM));
		uw2.addReplicate(1, Wash.UW2, CellType.C3H, CellCompartment.NE, getRemoteFile(s6_Uw2_NE));
		exps_uw2.add(uw2);

		final Experiment uw3 = new Experiment("Uw3", Wash.UW3, CellType.C3H);
		uw3.addReplicate(1, Wash.UW3, CellType.C3H, CellCompartment.CM, getRemoteFile(s12_Uw3_CM));
		uw3.addReplicate(1, Wash.UW3, CellType.C3H, CellCompartment.NE, getRemoteFile(s12_Uw3_NE));
		exps_uw3.add(uw3);

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
					if (wash == Wash.PREW1) {
						System.out.println("asdf");
					}
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
		final List<Experiment> list = new ArrayList<Experiment>();
		list.addAll(exps_uw1);
		list.addAll(exps_uw2);
		list.addAll(exps_uw3);
		list.addAll(exps_cw1);
		list.addAll(exps_prew1);
		list.addAll(exps_prew2);
		return list;
	}

}
