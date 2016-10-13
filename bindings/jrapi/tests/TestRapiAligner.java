/************************************************************************************
 * This code is published under the The MIT License.
 *
 * Copyright (c) 2016 Center for Advanced Studies,
 *                      Research and Development in Sardinia (CRS4), Pula, Italy.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ************************************************************************************/



import it.crs4.rapi.*;
import it.crs4.rapi.RapiUtils;
import it.crs4.rapi.AlignOp;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.junit.*;
import static org.junit.Assert.*;

public class TestRapiAligner
{
  private static File jrapiBaseDir;

  private Opts rapiOpts;
  private Ref refObj;
  private Batch reads;
  private AlignerState aligner;
  private File miniRefPath;

  @BeforeClass
  public static void initSharedObj()
  {
    RapiUtils.loadPlugin();
  }

  @Before
  public void init() throws RapiException, IOException
  {
    rapiOpts = new Opts();
    rapiOpts.setShareRefMem(false); // for travis
    Rapi.init(rapiOpts);

    reads = new Batch(2);
    TestUtils.appendSeqsToBatch(TestUtils.readMiniRefSeqs(), reads);

    // load the reference
    miniRefPath = new File(jrapiBaseDir, TestUtils.RELATIVE_MINI_REF);
    refObj = new Ref(miniRefPath.getAbsolutePath());

    aligner = new AlignerState(rapiOpts);
    aligner.alignReads(refObj, reads);
  }

  @After
  public void tearDown() throws RapiException
  {
    aligner = null;
    refObj.unload();
    refObj = null;
    if (reads != null)
      reads.clear();
    reads = null;
    Rapi.shutdown();
  }

  @Test
  public void testReadAttributes() throws RapiException
  {
    Read rapiRead = reads.getRead(0, 0);
    assertEquals("read_00", rapiRead.getId());
    // shortcut attributes that access the first Alignment
    assertFalse(rapiRead.getPropPaired());
    assertTrue(rapiRead.getMapped());
    assertFalse(rapiRead.getReverseStrand());
    assertEquals(60, rapiRead.getMapq());
    assertEquals(60, rapiRead.getScore());
  }

  @Test
  public void testAlignmentStruct() throws RapiException
  {
    Read rapiRead = reads.getRead(0, 0);
    assertEquals("read_00", rapiRead.getId());
    assertTrue(rapiRead.getNAlignments() > 0);

    Alignment aln = rapiRead.getAln(0);
    assertEquals("chr1", aln.getContig().getName());
    assertEquals(32461, aln.getPos());
    assertEquals("60M", aln.getCigarString());

    AlignOp[] ops = aln.getCigarOps();
    assertEquals(1, ops.length);
    assertEquals(AlignOp.Type.Match, ops[0].getType());
    assertEquals(60, ops[0].getLen());

    assertEquals(60, aln.getMapq());
    assertEquals(60, aln.getScore());

    // Flag 65 is 'p1'
    assertTrue(aln.getPaired());
    assertFalse(aln.getPropPaired());
    assertTrue(aln.getMapped());
    assertFalse(aln.getReverseStrand());
    assertFalse(aln.getSecondaryAln());
    assertEquals(0, aln.getNMismatches());
    assertEquals(0, aln.getNGapOpens());
    assertEquals(0, aln.getNGapExtensions());

    HashMap<String,Object> expected_tags = new HashMap<String, Object>();
    expected_tags.put("MD", "60");
    expected_tags.put("XS", 0L);
    // the SAM has more tags (NM and AS), but RAPI makes that information
    // available through other members of the Alignment structure.

    assertEquals(expected_tags, aln.getTags());

    // check out reads that don't align perfectly
    aln = reads.getRead(1, 0).getAln(0);
    assertEquals(60, aln.getMapq());
    assertEquals(51, aln.getScore());
    assertEquals("11M3D49M", aln.getCigarString());

    assertArrayEquals(
        new AlignOp[] { new AlignOp('M', 11), new AlignOp('D', 3), new AlignOp('M', 49) },
        aln.getCigarOps());

    assertEquals(3, aln.getNMismatches());

    String md_tag = (String)aln.getTags().get("MD");
    assertEquals("11^CCC49", md_tag);

    assertEquals(0, aln.getNGapOpens());
    assertEquals(0, aln.getNGapExtensions());

    aln = reads.getRead(2, 0).getAln(0);
    assertEquals(60, aln.getMapq());
    assertEquals(48, aln.getScore());
    assertEquals("13M3I44M", aln.getCigarString());

    assertArrayEquals(
        new AlignOp[] { new AlignOp('M', 13), new AlignOp('I', 3), new AlignOp('M', 44) },
        aln.getCigarOps());
    assertEquals(3, aln.getNMismatches());

    md_tag = (String)aln.getTags().get("MD");
    assertEquals("57", md_tag);

    aln = reads.getRead(3, 0).getAln(0);
    assertEquals(60, aln.getMapq());
    assertEquals(50, aln.getScore());
    assertEquals("60M", aln.getCigarString());
    assertEquals(2, aln.getNMismatches());
    md_tag = (String)aln.getTags().get("MD");
    assertEquals("15T16C27", md_tag);
  }

  @Test
  public void testGetInsertSize()
  {
    Alignment aln1 = reads.getRead(0, 0).getAln(0);
    Alignment aln2 = reads.getRead(0, 1).getAln(0);

    assertEquals(121L, Rapi.getInsertSize(aln1, aln2));
    assertEquals(-121L, Rapi.getInsertSize(aln2, aln1));
    assertEquals(0L, Rapi.getInsertSize(aln1, aln1));
  }

  @Test
  public void testAlnIterator() throws RapiException
  {
    Read rapiRead = reads.getRead(0, 0);

    Iterator<Alignment> it = rapiRead.iterator();
    assertTrue(it.hasNext());

    Alignment aln = it.next();
    assertEquals("chr1", aln.getContig().getName());
    assertEquals(32461, aln.getPos());
    assertEquals("60M", aln.getCigarString());

    assertFalse(it.hasNext());
  }

  public static void main(String args[])
  {
    TestUtils.testCaseMainMethod(TestRapiAligner.class.getName(), args);
  }
}

// vim: set et sw=2
