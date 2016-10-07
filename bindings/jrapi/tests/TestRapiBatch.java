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

import org.junit.*;
import static org.junit.Assert.*;

import java.util.Iterator;
import java.util.List;

public class TestRapiBatch
{
  private Batch b;
  private List<String[]> someReads;

  @BeforeClass
  public static void initSharedObj()
  {
    RapiUtils.loadPlugin();
  }

  @Before
  public void init() throws RapiException
  {
    someReads = TestUtils.getSomeReads();

    Rapi.init(new Opts());
    b = new Batch(2);
  }

  @After
  public void tearDown() throws RapiException
  {
    Rapi.shutdown();
  }

  @Test
  public void testInit()
  {
    assertEquals(2, b.getNReadsPerFrag());
    assertEquals(0, b.getCapacity());
    assertEquals(0, b.getLength());
    assertEquals(0, b.getNFragments());
  }

  @Test
  public void testReserve() throws RapiException
  {
    assertEquals(0, b.getCapacity());
    b.reserve(4);
    assertEquals(2, b.getNReadsPerFrag());
    assertEquals(4, b.getCapacity());
    assertEquals(0, b.getLength());
    assertEquals(0, b.getNFragments());
  }

  @Test(expected=RapiOutOfMemoryError.class)
  public void testImpossibleReserve() throws RapiException
  {
    assertEquals(0, b.getCapacity());
    b.reserve(99999999999L);
  }


  @Test
  public void testAppend() throws RapiException
  {
    String[] aRead = someReads.get(0);

    assertEquals(2, b.getNReadsPerFrag());

    assertEquals(0, b.getLength());
    assertEquals(0, b.getNFragments());

    b.append(aRead[0], aRead[1], aRead[2], Rapi.QENC_SANGER);

    assertEquals(1, b.getLength());
    assertEquals(0, b.getNFragments());

    b.append(aRead[0], aRead[3], aRead[4], Rapi.QENC_SANGER);

    assertEquals(2, b.getLength());
    assertEquals(1, b.getNFragments());

    assertTrue(2 <= b.getCapacity());
  }

  @Test(expected=RapiOutOfMemoryError.class)
  public void testImpossibleAllocation() throws RapiException
  {
    b.reserve(2000000000);
  }

  @Test
  public void testGetRead() throws RapiException
  {
    String[] aRead = someReads.get(0);
    b.append(aRead[0], aRead[1], aRead[2], Rapi.QENC_SANGER);
    b.append(aRead[0], aRead[3], aRead[4], Rapi.QENC_SANGER);

    Read r_get = b.getRead(0, 0);

    assertNotNull(r_get);
    assertEquals(aRead[0], r_get.getId());
    assertEquals(aRead[1], r_get.getSeq());
    assertEquals(aRead[2], r_get.getQual());
    assertEquals(aRead[1].length(), r_get.getLength());

    r_get = b.getRead(0, 1);
    assertNotNull(r_get);
    assertEquals(aRead[0], r_get.getId());
    assertEquals(aRead[3], r_get.getSeq());
    assertEquals(aRead[4], r_get.getQual());
  }

  @Test(expected=RapiInvalidParamException.class)
  public void testGetError1() throws RapiException
  {
    // verify out-of-bounds checking
    loadSomeReads(1);
    Read r = b.getRead(1, 0);
  }

  @Test(expected=RapiInvalidParamException.class)
  public void testGetError2() throws RapiException
  {
    // verify out-of-bounds checking
    loadSomeReads(1);
    Read r = b.getRead(-1, 0);
  }

  @Test(expected=RapiInvalidParamException.class)
  public void testGetError3() throws RapiException
  {
    // verify out-of-bounds checking
    loadSomeReads(1);
    Read r = b.getRead(0, 2);
  }

  @Test(expected=RapiInvalidParamException.class)
  public void testGetError4() throws RapiException
  {
    // verify out-of-bounds checking
    loadSomeReads(1);
    Read r = b.getRead(0, -1);
  }

  @Test
  public void testFormatSAMBatch() throws RapiException
  {
    loadSomeReads(2);
    String output = Rapi.formatSamBatch(b);
    assertTrue(output.length() > 0);

    String[] lines = output.split("\n");
    assertEquals(4, lines.length);
    // compare Read ids.  Each element in someReads contains two reads, so two lines of SAM
    for (int i = 0; i < lines.length / 2; ++i) {
      assertEquals(someReads.get(i)[0], lines[2*i].split("\t")[0]);
      assertEquals(someReads.get(i)[0], lines[2*i+1].split("\t")[0]);
    }
  }

  @Test
  public void testFormatSAMBatchIndexed() throws RapiException
  {
    loadSomeReads(2);
    String output = Rapi.formatSamBatch(b, 1); // generate SAM for the 2nd fragment in the batch
    assertTrue(output.length() > 0);

    String[] lines = output.split("\n");
    assertEquals(2, lines.length);
    // compare Read ids
    for (int i = 0; i < lines.length; ++i) {
      assertEquals(someReads.get(1)[0], lines[i].split("\t")[0]);
    }
  }

  @Test
  public void testIterator() throws RapiException
  {
    loadSomeReads(2);
    String[] firstReads = someReads.get(0);

    Iterator<Fragment> it = b.iterator();
    assertTrue(it.hasNext());

    Fragment f = it.next();
    assertEquals(b.getNReadsPerFrag(), f.getLength());

    assertEquals(firstReads[0], f.get(0).getId());
    assertEquals(firstReads[0], f.get(1).getId());
    assertEquals(firstReads[1], f.get(0).getSeq());
    assertEquals(firstReads[2], f.get(0).getQual());
    assertEquals(firstReads[3], f.get(1).getSeq());
    assertEquals(firstReads[4], f.get(1).getQual());

    assertTrue(it.hasNext());
    it.next();
    assertFalse(it.hasNext());
  }

  @Test
  public void testFragment() throws RapiException
  {
    loadSomeReads(1);
    String[] firstReads = someReads.get(0);

    Fragment f = b.iterator().next();
    Iterator<Read> fit = f.iterator();

    Read read1 = fit.next();
    assertEquals(firstReads[0], read1.getId());
    assertEquals(firstReads[1], read1.getSeq());
    assertEquals(firstReads[2], read1.getQual());

    Read read2 = fit.next();
    assertEquals(firstReads[0], read2.getId());
    assertEquals(firstReads[3], read2.getSeq());
    assertEquals(firstReads[4], read2.getQual());

    assertFalse(fit.hasNext());
  }

  private void loadSomeReads(int n_fragments) throws RapiException
  {
    if (n_fragments > someReads.size())
      throw new IllegalArgumentException("Requested n_fragments (" + n_fragments + ") is greated than number of reads available for test");

    for (int i = 0; i < n_fragments; ++i) {
      String[] fragment = someReads.get(i);
      b.append(fragment[0], fragment[1], fragment[2], Rapi.QENC_SANGER);
      b.append(fragment[0], fragment[3], fragment[4], Rapi.QENC_SANGER);
    }
  }

  public static void main(String args[])
  {
    TestUtils.testCaseMainMethod(TestRapiBatch.class.getName(), args);
  }
}

// vim: set et sw=2
