package edu.scripps.yates.nucleome.turboID.string;

import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.hash.THashSet;

/**
 * Reads the string interactor files from our baits and:<br>
 * generate the input file with weights 1.<br>
 * generates the functional annotation file.<br>
 * 
 * @author salvador
 *
 */
public class SafeInputFileGenerator {
	private File[] stringFiles = null;
	private final File inputFolder;
	private final UniprotProteinLocalRetriever uplr;
	private final boolean useStringScoreAsWeight;

	/**
	 * 
	 * @param args: full path to folder in which string files are<br>
	 *              full path to uniprot folder
	 */
	public static void main(String[] args) {
		final File folder = new File(args[0]);
		final File[] stringFiles = folder.listFiles(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String name) {
				if (name.startsWith("string_interactions")) {
					return true;
				}
				return false;
			}
		});
		final File uniprotFolder = new File(args[1]);
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotFolder, true);

		boolean useStringScoreAsWeight = false;
		if (args.length > 2) {
			useStringScoreAsWeight = Boolean.valueOf(args[2]);
		}
		final SafeInputFileGenerator sifg = new SafeInputFileGenerator(stringFiles, uplr, useStringScoreAsWeight);
		try {
			sifg.run();
			System.out.println("Everything ok!!");
		} catch (final IOException e) {
			e.printStackTrace();
		}
	}

	public SafeInputFileGenerator(File[] stringFiles, UniprotProteinLocalRetriever uplr,
			boolean useStringScoreAsWeight) {
		this.stringFiles = stringFiles;
		this.inputFolder = stringFiles[0].getParentFile();
		this.uplr = uplr;
		this.useStringScoreAsWeight = useStringScoreAsWeight;
	}

	private File getNetworkFile() {
		return new File(this.inputFolder.getAbsolutePath() + File.separator + "safe_input_network_file.tsv");
	}

	private File getFunctionalAnnotationFile() {
		return new File(
				this.inputFolder.getAbsolutePath() + File.separator + "safe_input_functional_annotation_file.tsv");
	}

	public void run() throws IOException {
		createNetworkFile(getNetworkFile());
	}

	private void createNetworkFile(File networkFile) throws IOException {
		FileWriter fw = null;
		try {
			fw = new FileWriter(networkFile);
			final Set<String> interactions = new THashSet<String>();
			int redundantInteractions = 0;

			for (final File stringFile : stringFiles) {
				final List<String> lines = Files.readAllLines(stringFile.toPath());
				final TObjectIntMap<String> indexByHeader = getIndexByHeader(lines.get(0));
				for (int row = 1; row < lines.size(); row++) {
					final String line = lines.get(row);
					final String[] split = line.split("\t");
					final String gene1 = split[indexByHeader.get("#node1")];
					final String gene2 = split[indexByHeader.get("node2")];
					final String interactionKey = getInteractionKey(gene1, gene2);
					if (interactions.contains(interactionKey)) {
						redundantInteractions++;
						continue;
					} else {
						interactions.add(interactionKey);
					}
					double weigth = 1;
					if (useStringScoreAsWeight) {
						weigth = Double.valueOf(split[indexByHeader.get("combined_score")]);
					}
					fw.write(gene1 + "\t" + gene2 + "\t" + weigth + "\n");
				}
			}
			System.out.println(interactions.size() + " unique PPIs from " + stringFiles.length + " input files ("
					+ redundantInteractions + " were redundant)");
		} finally {
			fw.close();
			System.out.println("File created at '" + networkFile.getAbsolutePath() + "'");
		}
	}

	private String getInteractionKey(String gene1, String gene2) {
		if (gene1.compareTo(gene2) > 0) {
			return gene1 + gene2;
		} else {
			return gene2 + gene1;
		}
	}

	private TObjectIntMap<String> getIndexByHeader(String firstLine) {
		final TObjectIntMap<String> ret = new TObjectIntHashMap<String>();
		final String[] split = firstLine.split("\t");
		int index = 0;
		for (final String header : split) {
			ret.put(header, index++);
		}
		return ret;
	}
}
