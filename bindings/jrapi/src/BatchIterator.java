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



package it.crs4.rapi;

public class BatchIterator implements java.util.Iterator<Fragment> {
  private Batch batch;
  private long nextFragment;
  private long end;
  private Fragment frag;

  public BatchIterator(Batch b)
  {
    if (b == null)
      throw new NullPointerException();
    batch = b;
    nextFragment = 0;
    end = batch.getNFragments();
    frag = new Fragment(batch.getNReadsPerFrag());
  }

  public Fragment next()
  {
    if (nextFragment >= end)
      throw new java.util.NoSuchElementException();

    for (int i = 0; i < batch.getNReadsPerFrag(); ++i)
      frag.set(i, batch.getRead(nextFragment, i));

    nextFragment += 1;
    return frag;
  }

  public boolean hasNext()
  {
    return nextFragment < end;
  }

  public void remove()
  {
    throw new UnsupportedOperationException();
  }
}
