package edu.scripps.yates.nucleome.turboID;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.set.hash.THashSet;

/**
 * This class has method to generate the required table to use in as input with
 * UMAP clustering
 * 
 * @author salvador
 *
 */
public class UMAPClustering {

	private final File folder;
	private final List<TurboIDExperiment> experiments = new ArrayList<TurboIDExperiment>();
	private File outputFile;

	private final static Logger log = Logger.getLogger(UMAPClustering.class);
	//////////////////////////
	//
	private final boolean useDistributedIntensity = true;
	private final boolean useNormalizedIntensity = false;

	//
	///////////////////////
	public UMAPClustering(File folder, TurboIDExperiment cyExperiment, TurboIDExperiment nuExperiment)
			throws IOException {
		this.folder = folder;
		if (nuExperiment == null || cyExperiment == null) {
			throw new IllegalArgumentException("control and real experiment need to be provided");
		}
		this.experiments.add(nuExperiment); // first the nu
		this.experiments.add(cyExperiment); // sedond the cy
		if (!folder.exists()) {
			folder.mkdirs();
		}

	}

	public void run(boolean onlyTM, boolean onlyNonTM) throws IOException {

		String infix = "norm";
		if (useDistributedIntensity) {
			infix = "distr";
		}
		outputFile = new File(folder.getAbsoluteFile() + File.separator + "umap_" + infix + ".tsv");

		writeUMAPInputFile(onlyTM, onlyNonTM);

	}

	private String getBaitName(String name) {
		if (name.equals(TurboIDExperimentType.EMD.name())) {
			return "Emd";
		} else if (name.equals(TurboIDExperimentType.LBR.name())) {
			return "Lbr";
		} else if (name.equals(TurboIDExperimentType.MAN1.name())) {
			return "Man1";
		} else if (name.equals(TurboIDExperimentType.SUN1.name())) {
			return "Sun1";
		}
		return null;
	}

	private void writeUMAPInputFile(boolean onlyTM, boolean onlyNonTM) throws IOException {
		FileWriter fw = null;

		try {
			final File benchmarksFile = new File(
					folder.getParentFile().getAbsolutePath() + File.separator + "benchmarks_20191104.txt");
			final Map<String, String> benchmarks = Benchmarks.getBenchmarks(benchmarksFile);
			fw = new FileWriter(outputFile);

			final Set<String> accs = new THashSet<String>();
			for (final TurboIDExperiment ex : experiments) {
				accs.addAll(ex.keySet());
			}
			final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(
					TurboIDDataAnalysis.uniprotReleasesFolder, true);
			final Map<String, Entry> annotatedProteins = uplr.getAnnotatedProteins(null, accs);
			// header
			final List<TurboIDExperimentType> baits = TurboIDExperimentType.getBaits();
			for (final TurboIDExperimentType bait : baits) {
				fw.write("\t" + getBaitName(bait.name()) + "_Nu");
			}
			for (final TurboIDExperimentType bait : baits) {
				fw.write("\t" + getBaitName(bait.name()) + "_Cy");
			}
			fw.write("\tTM\tbenchmark\tlocalization");
			fw.write("\n");

			// data
			for (final String acc : accs) {
				boolean isvalid = false;
				boolean isTM = false;
				String gene = null;
				final List<Pair<String, String>> geneName = UniprotEntryUtil.getGeneName(annotatedProteins.get(acc),
						true, true);
				if (geneName != null && !geneName.isEmpty()) {
					gene = geneName.get(0).getFirstelement();
				}
				if (gene == null) {
					continue;
				}
				fw.write(gene);
				for (final TurboIDExperiment ex : experiments) {
					final ProteinFromTurboID protein = ex.get(acc);
					if (protein != null) {
						if (onlyTM && !protein.isTransmembrane()) {
							continue;
						}
						if (onlyNonTM && protein.isTransmembrane()) {
							continue;
						}
						isTM = protein.isTransmembrane();
					}
					isvalid = true;

					// Nu
					for (final TurboIDExperimentType bait : baits) {
						final TDoubleList toAverage = new TDoubleArrayList();
						for (final Replicate replicate : Replicate.values(ex.getFraction())) {
							// do not use Arerun
							if (replicate == Replicate.Arerun) {
								continue;
							}
							for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
								if (channel.getExpType() != bait) {
									continue;
								}

								if (channel.getReplicate() == replicate) {

									double intensity = Double.NaN;
									if (protein != null) {
										if (useDistributedIntensity) {
											intensity = protein.getDistributedIntensitiesWithNormalizedIntensities()
													.get(channel);
										}
										if (useNormalizedIntensity) {
											intensity = protein.getNormalizedIntensities(bait).get(channel);
										}
									}
									if (!Double.isNaN(intensity) && Double.compare(0.0, intensity) != 0) {
										toAverage.add(intensity);
									}

								}
							}
						}
						if (toAverage.size() > 0) {
							fw.write("\t" + Maths.mean(toAverage));
						} else {
							fw.write("\t" + 0.0);
						}
					}
				}
				if (isvalid) {
					fw.write("\t" + isTM);
					final boolean isBenchmark = benchmarks.containsKey(gene);
					fw.write("\t" + isBenchmark);
					if (isBenchmark) {
						fw.write("\t" + benchmarks.get(gene));
					} else {
						fw.write("\t");
					}
					fw.write("\n");
				}
			}

		} finally {
			if (fw != null) {
				fw.close();
			}
		}
	}

}
