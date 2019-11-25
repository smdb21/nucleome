package edu.scripps.yates.nucleome.turboID.comparison;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.math.NumberUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.excel.ExcelReader;
import edu.scripps.yates.nucleome.turboID.TurboIDExperimentType;
import edu.scripps.yates.utilities.venndata.VennData;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.hash.THashSet;

/**
 * Having 2 different tables with the result of the TurboID ANOVA analysis, we
 * want to compare the data.
 * 
 * @author salvador
 *
 */
public class TurboIDComparison {
	private static final Logger log = Logger.getLogger(TurboIDComparison.class);
	// not optimized experiment
	private final File experiment1File = new File(
			"Z:\\share\\Salva\\data\\4D_Nucleome\\TurboID\\input\\NU_anova_analysis_filtered_DATA.xlsx");
	// optimized params experiment:
	private final File experiment2File = new File(
			"Z:\\share\\Salva\\data\\4D_Nucleome\\TurboID\\input\\optimized_params\\NuCy_anova_distribInt_optimized_params_DATA.xlsx");
	private final double pValueThreshold;
	private TIntIntHashMap correspondingClusters;
	private final boolean useCorrespondingClusters = false;
	private final boolean useMinValueClusterColumn = true;
	private final int minSPC = 3;

	public static void main(String[] args) {
		try {
			final double pvalueThreshold = 0.001;
			new TurboIDComparison(pvalueThreshold).run();
			System.exit(0);
		} catch (final IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	public TurboIDComparison(double pvalueThreshold) {
		pValueThreshold = pvalueThreshold;
	}

	private void run() throws IOException {
		FileWriter fw = null;
		FileWriter fw2 = null;
		try {
			fw = new FileWriter(new File(experiment2File.getParentFile().getAbsolutePath() + File.separator
					+ "Old_vs_Optimized_params.txt"));
			fw2 = new FileWriter(new File(experiment2File.getParentFile().getAbsolutePath() + File.separator
					+ "Old_vs_Optimized_params_correlation.txt"));
			fw.write("\tOLD\tOptimized\tintersection\n");

			final Experiment exp1 = readDataFromFile(experiment1File);
			final Experiment exp2 = readDataFromFile(experiment2File);
			final VennData vennProteinsUnderThreshold = new VennData("Proteins p<" + pValueThreshold, "OLD",
					exp1.getProteinsByAcc().keySet(), "Optimized", exp2.getProteinsByAcc().keySet(), null, null);
			final VennData vennProteinsClustered = new VennData("Proteins clustered", "OLD",
					exp1.getClusteredProteins(), "Optimized", exp2.getClusteredProteins(), null, null);
			System.out.println(vennProteinsClustered.getIntersectionsText("Proteins clustered"));
			System.out.println(vennProteinsUnderThreshold.getIntersectionsText("Proteins p<" + pValueThreshold));
			final VennData vennTotalProteins = new VennData("Total proteins", "OLD", exp1.getTotalProteins(),
					"Optimized", exp2.getTotalProteins(), null, null);
			fw.write("Num proteins p<" + pValueThreshold + "\t" + exp1.getNumSize() + "\t" + exp2.getNumSize() + "\t"
					+ vennProteinsUnderThreshold.getIntersection12().size() + "\n");
			fw.write("Num proteins clustered" + "\t" + exp1.getClusteredProteins().size() + "\t"
					+ exp2.getClusteredProteins().size() + "\t" + vennProteinsClustered.getIntersection12().size()
					+ "\n");
			System.out.println(vennTotalProteins.getIntersectionsText("Total proteins"));
			fw.write("Num Total Proteins\t" + exp1.getTotalProteins().size() + "\t" + exp2.getTotalProteins().size()
					+ "\t" + vennTotalProteins.getIntersection12().size() + "\n");
			fw.write("Num PSMs\t" + exp1.getTotalSPC() + "\t" + exp2.getTotalSPC() + "\n");
			fw.write("Num Complete Proteins\t" + exp1.getNumCompleteProteins() + "\t" + exp2.getNumCompleteProteins()
					+ "\n");

			fw.write("\n\n\n");

			final boolean[] truefalseArray = { true, false };
			for (final boolean useEuclideanDistance : truefalseArray) {
				if (useEuclideanDistance) {
					fw.write("Clusters between experiments mapped by Euclidean distance\n");
				} else {
					fw.write("Clusters between experiments mapped by Pearson Correlation\n");
				}
				fw.write("Cluster in exp1" + "\t" + "Cluster in exp2" + "\t" + "Num proteins in exp1" + "\t"
						+ "Num proteins in exp2" + "\t" + "intersection" + "\t");
				for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
					fw.write("Exp1_Avg_" + bait.name() + "\t");
					fw.write("Exp2_Avg_" + bait.name() + "\t");
				}
				fw.write("\n");

				for (int cluster = 1; cluster <= Math.max(exp1.getNumClusters(), exp2.getNumClusters()); cluster++) {
					fw.write(cluster + "\t");
					final TObjectDoubleMap<TurboIDExperimentType> averageQuantsByCluster = exp1
							.getAverageQuantsByCluster(cluster);
					int cluster2 = getCorrespondingCluster(cluster);
					if (cluster2 < 1) {
						if (useEuclideanDistance) {
							cluster2 = exp2.getCloserClusterByEuclideanDistance(averageQuantsByCluster);
						} else {
							cluster2 = exp2.getCloserClusterByPearsonCorrelation(averageQuantsByCluster);
						}
					}

					fw.write(cluster2 + "\t");
					fw.write(exp1.getProteinsByClusterNumber(cluster).size() + "\t"
							+ exp2.getProteinsByClusterNumber(cluster2).size() + "\t");
					final Set<String> proteins1 = exp1.getProteinsByClusterNumber(cluster).stream().map(p -> p.getAcc())
							.collect(Collectors.toSet());
					final Set<String> proteins2 = exp2.getProteinsByClusterNumber(cluster2).stream()
							.map(p -> p.getAcc()).collect(Collectors.toSet());
					final VennData venn2 = new VennData(cluster + " vs " + cluster2, "" + cluster, proteins1,
							"" + cluster2, proteins2, null, null);
					final int intersection = venn2.getIntersection12().size();
					fw.write(intersection + "\t");

					final TObjectDoubleMap<TurboIDExperimentType> averageQuantsByCluster2 = exp2
							.getAverageQuantsByCluster(cluster2);
					for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
						final double avg1 = averageQuantsByCluster.get(bait);
						final double avg2 = averageQuantsByCluster2.get(bait);
						fw.write(avg1 + "\t" + avg2 + "\t");
					}
					fw.write("\n");
				}

				fw2.write("Protein\t");
				for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
					fw2.write("Exp1_" + bait.name() + "\t");
					fw2.write("Exp2_" + bait.name() + "\t");
				}
				for (final ProteinToCompare protein : exp1.getProteins()) {
					if (exp2.getProteinsByAcc().containsKey(protein.getAcc())) {
						fw2.write(protein.getAcc() + "\t");
						final ProteinToCompare protein2 = exp2.getProteinsByAcc().get(protein.getAcc());
						for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
							double quant = protein.getQuant(bait);
							if (Double.isNaN(quant)) {
								quant = 0.0;
							}
							double quant2 = protein2.getQuant(bait);
							if (Double.isNaN(quant2)) {
								quant2 = 0.0;
							}
							fw2.write(quant + "\t" + quant2 + "\t");
						}
						fw2.write("\n");
					}
				}
			}
		} finally {
			fw.close();
			fw.close();
		}
	}

	private int getCorrespondingCluster(int cluster) {
		if (correspondingClusters == null) {
			correspondingClusters = new TIntIntHashMap();
			if (useCorrespondingClusters) {
				correspondingClusters.put(1, 1);
				correspondingClusters.put(2, 2);
				correspondingClusters.put(3, 10);
				correspondingClusters.put(4, 5);
				correspondingClusters.put(5, 11);
				correspondingClusters.put(6, 3);
				correspondingClusters.put(7, 8);
				correspondingClusters.put(8, 9);
				correspondingClusters.put(9, 4);
				correspondingClusters.put(10, 7);
				correspondingClusters.put(11, 6);
			}
		}
		return correspondingClusters.get(cluster);

	}

	private Experiment readDataFromFile(File file) throws IOException {
		final Experiment experiment = new Experiment();
		final ExcelReader excel1 = new ExcelReader(file, 0, 0);
		final int emdColumnIndex = excel1.getColumnIndex(0, "distNormInt EMD");
		final int lbrColumnIndex = excel1.getColumnIndex(0, "distNormInt LBR");
		final int sun1ColumnIndex = excel1.getColumnIndex(0, "distNormInt SUN1");
		final int man1ColumnIndex = excel1.getColumnIndex(0, "distNormInt MAN1");
		final int accColumnIndex = excel1.getColumnIndex(0, "accession");
		final int geneColumnIndex = excel1.getColumnIndex(0, "gene");
		final int rep1SPC = excel1.getColumnIndex(0, "SPC_A");
		final int rep2SPC = excel1.getColumnIndex(0, "SPC_B");
		final int rep3SPC = excel1.getColumnIndex(0, "SPC_Brerun");
		final int descriptionColumnIndex = excel1.getColumnIndex(0, "description");
		String columnName = "cluster p<" + pValueThreshold;
		if (useMinValueClusterColumn) {
			columnName += " & min_value";
		}
		final int clusterColumnIndex = excel1.getColumnIndex(0, columnName);
		final int pvalueColumnIndex = excel1.getColumnIndex(0, "distNormInt ANOVA p-value ");
		final int validColumnIndex = excel1.getColumnIndex(0, "is complete");
		final int totalSPCColumnIndex = excel1.getColumnIndex(0, "total SPCs");
		int row = 1;
		int totalSPCInExp = 0;
		int numValid = 0;
		final Set<String> clusteredProteins = new THashSet<String>();
		final Set<String> proteins = new THashSet<String>();
		while (true) {
			try {
				final String acc = excel1.getStringValue(0, row, accColumnIndex);
				if (acc == null) {
					log.info("No more data at row " + row);
					break;
				}
				proteins.add(acc);
				final String valid = excel1.getStringValue(0, row, validColumnIndex);
				if (!valid.equals("true")) {
					continue;
				}
				numValid++;
				final String gene = excel1.getStringValue(0, row, geneColumnIndex);
				final String description = excel1.getStringValue(0, row, descriptionColumnIndex);
				final double emd = Double.valueOf(excel1.getNumberValue(0, row, emdColumnIndex));
				final double lbr = Double.valueOf(excel1.getNumberValue(0, row, lbrColumnIndex));
				final double sun1 = Double.valueOf(excel1.getNumberValue(0, row, sun1ColumnIndex));
				final double man1 = Double.valueOf(excel1.getNumberValue(0, row, man1ColumnIndex));

				final double spcRep1 = Double.valueOf(excel1.getNumberValue(0, row, rep1SPC));
				final double spcRep2 = Double.valueOf(excel1.getNumberValue(0, row, rep2SPC));
				final double spcRep3 = Double.valueOf(excel1.getNumberValue(0, row, rep3SPC));
				if (spcRep1 < minSPC || spcRep2 < minSPC || spcRep3 < minSPC) {
					continue;
				}
				final double totalSPC = Double.valueOf(excel1.getNumberValue(0, row, totalSPCColumnIndex));
				totalSPCInExp += totalSPC;
				final double pvalue = Double.valueOf(excel1.getNumberValue(0, row, pvalueColumnIndex));
				if (Double.isNaN(pvalue) || pvalue > pValueThreshold) {
					continue;
				}
				final String numberValueString = excel1.getStringValue(0, row, clusterColumnIndex);
				Integer clusterNum = null;
				if (NumberUtils.isNumber(numberValueString)) {
					clusterNum = Double.valueOf(numberValueString).intValue();
					clusteredProteins.add(acc);
				} else if (pValueThreshold >= 0 && !useMinValueClusterColumn) {
					throw new IllegalArgumentException(
							"Error, all proteins with pvalue under threshold and complete should have a cluster number associated");
				}
				final ProteinToCompare protein = new ProteinToCompare(acc, gene, description, clusterNum, pvalue);
				protein.addQuant(TurboIDExperimentType.EMD, emd);
				protein.addQuant(TurboIDExperimentType.LBR, lbr);
				protein.addQuant(TurboIDExperimentType.MAN1, man1);
				protein.addQuant(TurboIDExperimentType.SUN1, sun1);
				experiment.addProtein(protein);
			} finally {
				row++;
			}
		}
		experiment.setTotalSPC(totalSPCInExp);
		experiment.setTotalProteins(proteins);
		experiment.setClusteredProteins(clusteredProteins);
		experiment.setNumCompleteProteins(numValid);
		log.info(experiment.getNumSize() + " proteins with p<" + pValueThreshold + ", " + proteins.size()
				+ " proteins in total and " + experiment.getTotalSPC() + " SPC in total "
				+ FilenameUtils.getBaseName(file.getAbsolutePath()));
		return experiment;
	}

}
