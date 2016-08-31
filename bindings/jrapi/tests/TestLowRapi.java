
import it.crs4.rapi.lowrapi.*;
import it.crs4.rapi.RapiUtils;

import java.io.File;

import org.junit.*;
import static org.junit.Assert.*;

public class TestLowRapi
{
  private static final String RELATIVE_MINI_REF = "../../tests/mini_ref/mini_ref.fasta";
  private static File jrapiBaseDir;

  private opts rapiOpts;
  private int error;
  private File miniRefPath;

  @BeforeClass
  public static void initSharedObj()
  {
    RapiUtils.loadPlugin(null);
  }

  @Before
  public void init()
  {
    miniRefPath = new File(jrapiBaseDir, RELATIVE_MINI_REF);
    rapiOpts = new opts();
    error = Rapi.init(rapiOpts);
    assertEquals(Rapi.NO_ERROR, error);
  }

  @After
  public void tearDown()
  {
    error = Rapi.shutdown();
    assertEquals(Rapi.NO_ERROR, error);
  }

  @Test
  public void testVersions()
  {
    assertTrue(Rapi.aligner_name().startsWith("bwa-mem"));
    assertTrue(Rapi.aligner_version().startsWith("0."));
    assertTrue(Rapi.plugin_version().startsWith("0."));
  }

  @Test
  public void testLoadUnloadRef()
  {
    ref refObj = new ref();
    error = Rapi.ref_load(miniRefPath.getAbsolutePath(), refObj);
    assertEquals(Rapi.NO_ERROR, error);

    assertEquals(1, refObj.getN_contigs());

    error = Rapi.ref_free(refObj);
    assertEquals(Rapi.NO_ERROR, error);
  }

  public static void main(String args[])
  {
    if (args.length == 1)
      RapiUtils.loadPlugin(args[0]);
    else
      RapiUtils.loadPlugin(null);

    org.junit.runner.JUnitCore.main(TestLowRapi.class.getName());
  }
}

// vim: set et sw=2
