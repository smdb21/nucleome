package edu.scripps.yates.nucleome.turboID;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.files.FileUtils;
import gnu.trove.map.hash.THashMap;

/**
 * This class reads into a file in which we have a table with the annotations
 * found by the Gingrass paper
 * 
 * @author salvador
 *
 */
public class Gingrass {
	private static final Logger log = Logger.getLogger(Gingrass.class);

	/**
	 * It returns a Map in which the key is the gene name and the value is the
	 * localization based on MMF
	 * 
	 * @param gingrassFile
	 * @return
	 * @throws IOException
	 */
	public static Map<String, String> getGingrassMMF(File gingrassFile) throws IOException {

		final List<String> genes = FileUtils.readColumnFromTextFile(gingrassFile, "\t", 1, 1);
		final List<String> mmf = FileUtils.readColumnFromTextFile(gingrassFile, "\t", 2, 1);
		final Map<String, String> ret = new THashMap<String, String>();
		for (int i = 0; i < mmf.size(); i++) {
			if (mmf.get(i) != null) {
				ret.put(genes.get(i), mmf.get(i));
			}
		}
		log.info(ret.size() + " proteins in Gingrass list (MMF localization) list from file "
				+ gingrassFile.getAbsolutePath());
		return ret;
	}

	/**
	 * It returns a Map in which the key is the gene name and the value is the
	 * localization based on MMF
	 * 
	 * @param gingrassFile
	 * @return
	 * @throws IOException
	 */
	public static Map<String, String> getGingrassSAFE(File gingrassFile) throws IOException {

		final List<String> genes = FileUtils.readColumnFromTextFile(gingrassFile, "\t", 1, 1);
		final List<String> safe = FileUtils.readColumnFromTextFile(gingrassFile, "\t", 3, 1);
		final Map<String, String> ret = new THashMap<String, String>();
		for (int i = 0; i < safe.size(); i++) {
			if (safe.get(i) != null) {
				ret.put(genes.get(i), safe.get(i));
			}
		}
		log.info(ret.size() + " proteins in Gingrass list (SAFE localization) from file "
				+ gingrassFile.getAbsolutePath());
		return ret;
	}
}
