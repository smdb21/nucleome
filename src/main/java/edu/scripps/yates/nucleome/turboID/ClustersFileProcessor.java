package edu.scripps.yates.nucleome.turboID;

import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FilenameUtils;

import com.google.common.io.Files;

import gnu.trove.map.hash.THashMap;

/**
 * This program parses the CSV file that is generated in R after printing the
 * cluster labels in a hierarchical cluster
 * 
 * @author salvador
 *
 */
public class ClustersFileProcessor {

	public static void main(String[] args) {
		final ClustersFileProcessor c = new ClustersFileProcessor();
		try {
			c.run("Z:\\share\\Salva\\data\\4D_Nucleome\\Xi data\\Xi data Feb 2019\\output");
		} catch (final IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	private void run(String folderPath) throws IOException {
		final File folder = new File(folderPath);
		final String[] list = folder.list(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String name) {
				if (name.contains("clusters.csv")) {
					return true;
				}
				return false;
			}
		});
		for (final String fileName : list) {
			final Map<Integer, List<String>> clusters = new THashMap<Integer, List<String>>();
			final File inputFile = new File(folder.getAbsolutePath() + File.separator + fileName);
			final List<String> readLines = Files.readLines(inputFile, Charset.defaultCharset());
			for (final String line : readLines) {
				final String[] split = line.split(",");
				final String gene = split[0].replace("\"", "");
				if (gene.equals("")) {
					continue;
				}
				if (split.length != 2) {
					System.out.println(split);
				}
				final int clusterID = Integer.valueOf(split[1].replace("\"", ""));
				if (clusters.containsKey(clusterID)) {
					clusters.get(clusterID).add(gene);
				} else {
					final List<String> list2 = new ArrayList<String>();
					list2.add(gene);
					clusters.put(clusterID, list2);
				}

			}
			final File outputFile = new File(
					folder.getAbsolutePath() + File.separator + FilenameUtils.getBaseName(fileName) + "_inColumns.tsv");
			final FileWriter fw = new FileWriter(outputFile);

			for (int i = 1; i < Integer.MAX_VALUE; i++) {
				if (!clusters.containsKey(i)) {
					break;
				}
				fw.write("Cluster " + i + "\t");
				final List<String> genes = clusters.get(i);
				Collections.sort(genes);
			}
			fw.write("\n");
			for (int index = 0; index < Integer.MAX_VALUE; index++) {
				boolean found = false;
				for (int clusterID = 1; clusterID < Integer.MAX_VALUE; clusterID++) {
					if (!clusters.containsKey(clusterID)) {
						break;
					}
					final List<String> genes = clusters.get(clusterID);
					if (genes.size() > index) {
						found = true;
						fw.write(genes.get(index));
					}
					fw.write("\t");
				}
				if (!found) {
					break;
				}
				fw.write("\n");
			}
			fw.close();
		}

	}
}
