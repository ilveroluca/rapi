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

import java.util.ArrayList;

public class Fragment implements Iterable<Read> {

  private ArrayList<Read> reads;

  public Fragment(int size) {
    reads = new ArrayList<Read>(size);
    for (int i = 0; i < size; ++i)
      reads.add(null);
  }

  public Read get(int i) {
    return reads.get(i);
  }

  public Read set(int i, Read r)
  {
    if (r == null)
      throw new IllegalArgumentException("Read cannot be null");

    return reads.set(i, r);
  }

  public java.util.Iterator<Read> iterator() {
    return reads.iterator();
  }

  public int getLength() {
    return reads.size();
  }
}
