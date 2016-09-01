
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
  private File miniRefPath;

  @BeforeClass
  public static void initSharedObj()
  {
    RapiUtils.loadPlugin(null);
  }

  @Before
  public void init() throws RapiException
  {
    miniRefPath = new File(jrapiBaseDir, RELATIVE_MINI_REF);
    rapiOpts = new opts();
    Rapi.init(rapiOpts);
  }

  @After
  public void tearDown() throws RapiException
  {
    Rapi.shutdown();
  }

  @Test
  public void testVersions()
  {
    assertTrue(Rapi.aligner_name().startsWith("bwa-mem"));
    assertTrue(Rapi.aligner_version().startsWith("0."));
    assertTrue(Rapi.plugin_version().startsWith("0."));
  }

  @Test
  public void testLoadUnloadRef() throws RapiException
  {
    ref refObj = new ref();
    Rapi.ref_load(miniRefPath.getAbsolutePath(), refObj);

    assertEquals(1, refObj.getN_contigs());

    Rapi.ref_free(refObj);
  }


  @Test
  public void testInstantiateAligner() throws RapiException
  {
    aligner_state state = Rapi.aligner_state_init(new opts());
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
