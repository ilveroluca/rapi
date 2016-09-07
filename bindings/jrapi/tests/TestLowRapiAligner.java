
import it.crs4.rapi.lowrapi.*;
import it.crs4.rapi.RapiUtils;
import it.crs4.rapi.AlignOp;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import org.junit.*;
import static org.junit.Assert.*;

public class TestLowRapiAligner
{
  private static File jrapiBaseDir;

  private opts rapiOpts;
  private ref refObj;
  private batch reads;
  private aligner_state aligner;
  private File miniRefPath;

  @BeforeClass
  public static void initSharedObj()
  {
    RapiUtils.loadPlugin(null);
  }

  @Before
  public void init() throws RapiException, IOException
  {
    rapiOpts = new opts();
    Rapi.init(rapiOpts);

    loadReads();

    // load the reference
    miniRefPath = new File(jrapiBaseDir, TestUtils.RELATIVE_MINI_REF);
    refObj = new ref();
    Rapi.ref_load(miniRefPath.getAbsolutePath(), refObj);

    aligner = Rapi.aligner_state_init(rapiOpts);
    Rapi.align_reads(refObj, reads, 0, reads.getN_frags(), aligner);
    System.err.println("============ Finished aligning reads");
  }

  @After
  public void tearDown() throws RapiException
  {
    // XXX: ************  need to delete aligner

    Rapi.ref_free(refObj);
    Rapi.reads_free(reads);
    Rapi.shutdown();
  }

  private void loadReads() throws RapiException, FileNotFoundException, IOException
  {
    List<String[]> someReads = TestUtils.readMiniRefSeqs();
    final int n_fragments = someReads.size();

    reads = new batch();
    Rapi.reads_alloc(reads, 2, n_fragments);

    for (int i = 0; i < n_fragments; ++i) {
      String[] pair = someReads.get(i);
      Rapi.set_read(reads, i, 0, pair[0], pair[1], pair[2], Rapi.QENC_SANGER);
      Rapi.set_read(reads, i, 1, pair[0], pair[3], pair[4], Rapi.QENC_SANGER);
    }
  }

  @Test
  public void testInstantiateAligner() throws RapiException
  {
    aligner_state state = Rapi.aligner_state_init(new opts());
  }

  @Test
  public void testReadAttributes() throws RapiException
  {
    read rapiRead = Rapi.get_read(reads, 0, 0);
    assertEquals("read_id:1", rapiRead.getId());
    // shortcut attributes that access the first alignment
    assertFalse(rapiRead.getPropPaired());
    assertTrue(rapiRead.getMapped());
    assertFalse(rapiRead.getReverseStrand());
    assertEquals(60, rapiRead.getMapq());
    assertEquals(60, rapiRead.getScore());
  }

  @Test
  public void testAlignmentStruct() throws RapiException
  {
    read rapiRead = Rapi.get_read(reads, 0, 0);
    assertEquals("read_id:1", rapiRead.getId());
    assertTrue(rapiRead.getNAlignments() > 0);

    alignment aln = rapiRead.getAln(0);
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
    assertEquals(0, aln.getN_mismatches());
    assertEquals(0, aln.getN_gap_opens());
    assertEquals(0, aln.getN_gap_extensions());

    /* Tags not yet implememnted
    expected_tags = dict(
        MD='60',
        XS=0)
# the SAM has more tags (NM and AS), but RAPI makes that information
# available through other members of the alignment structure.
      self.assertEqual(expected_tags, aln.get_tags())
      */

    // check out reads that don't align perfectly
    aln = Rapi.get_read(reads, 1, 0).getAln(0);
    assertEquals(60, aln.getMapq());
    assertEquals(51, aln.getScore());
    assertEquals("11M3D49M", aln.getCigarString());

    assertArrayEquals(
        new AlignOp[] { new AlignOp('M', 11), new AlignOp('D', 3), new AlignOp('M', 49) },
        aln.getCigarOps());

    assertEquals(3, aln.getN_mismatches());
    /* XXX: tags not implemented yet
    md_tag = aln.get_tags()['MD']
    assertEquals('11^CCC49', md_tag)
    */
    assertEquals(0, aln.getN_gap_opens());
    assertEquals(0, aln.getN_gap_extensions());

    aln = Rapi.get_read(reads, 2, 0).getAln(0);
    assertEquals(60, aln.getMapq());
    assertEquals(48, aln.getScore());
    assertEquals("13M3I44M", aln.getCigarString());

    assertArrayEquals(
        new AlignOp[] { new AlignOp('M', 13), new AlignOp('I', 3), new AlignOp('M', 44) },
        aln.getCigarOps());
    assertEquals(3, aln.getN_mismatches());

    /* XXX: tags not implemented yet
      md_tag = aln.get_tags()['MD']
      self.assertEqual('57', md_tag)
      */

    aln = Rapi.get_read(reads, 3, 0).getAln(0);
    assertEquals(60, aln.getMapq());
    assertEquals(50, aln.getScore());
    assertEquals("60M", aln.getCigarString());
    assertEquals(2, aln.getN_mismatches());
    /* XXX: tags not implemented yet
      md_tag = aln.get_tags()['MD'];
      assertEquals('15T16C27', md_tag);
      */
  }

  public static void main(String args[])
  {
    if (args.length == 1)
      RapiUtils.loadPlugin(args[0]);
    else
      RapiUtils.loadPlugin(null);

    org.junit.runner.JUnitCore.main(TestLowRapiAligner.class.getName());
  }
}

// vim: set et sw=2
