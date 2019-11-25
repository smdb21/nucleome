package edu.scripps.yates.nucleome.turboID;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.files.FileUtils;
import gnu.trove.map.hash.THashMap;

/**
 * This class reads into a file in which we have a table with the selected
 * protein benchmarks
 * 
 * @author salvador
 *
 */
public class Benchmarks {
	private static final Logger log = Logger.getLogger(Benchmarks.class);

	/**
	 * It returns a Map in which the key is the gene name and the value is the
	 * localization
	 * 
	 * @param benchmarksFile
	 * @return
	 * @throws IOException
	 */
	public static Map<String, String> getBenchmarks(File benchmarksFile) throws IOException {

		final List<String> genes = FileUtils.readColumnFromTextFile(benchmarksFile, "\t", 1, 2);
		final List<String> localization = FileUtils.readColumnFromTextFile(benchmarksFile, "\t", 5, 2);
		final Map<String, String> ret = new THashMap<String, String>();
		for (int i = 0; i < genes.size(); i++) {
			ret.put(genes.get(i), localization.get(i));
		}
		log.info(ret.size() + " proteins in benchmark list from file " + benchmarksFile.getAbsolutePath());
		return ret;
	}
}
