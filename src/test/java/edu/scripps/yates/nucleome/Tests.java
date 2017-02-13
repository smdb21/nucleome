package edu.scripps.yates.nucleome;

import org.junit.Test;

import edu.scripps.yates.nucleome.filters.GOFilter;

public class Tests {
	@Test
	public void testingGOFilter() {
		GOFilter filter = new GOFilter(" ");
		filter.isValid("O54962");
	}
}
