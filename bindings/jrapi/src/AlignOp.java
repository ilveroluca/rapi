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
import java.util.List;
import java.util.regex.*;

public class AlignOp {

  public static enum Type {
    Match,
    Insert,
    Delete,
    SoftClip,
    HardClip,
    Skip,
    Pad;

    public char getSymbol()
    {
      switch (this) {
        case Match:
          return 'M';
        case Insert:
          return 'I';
        case Delete:
          return 'D';
        case SoftClip:
          return 'S';
        case HardClip:
          return 'H';
        case Skip:
          return 'N';
        case Pad:
          return 'P';
        default:
          throw new RuntimeException("BUG!  Missing cigar operator in AlignOp::Type::getSymbol.");
      }
    }

    public static Type fromSymbol(String sym) {
      if (sym.length() != 1)
        throw new IllegalArgumentException("Unrecognized alignment operation symbol " + sym);
      else
        return fromSymbol(sym.charAt(0));
    }

    public static Type fromSymbol(char sym) {
      switch (sym) {
        case 'M':
          return Type.Match;
        case 'I':
          return Type.Insert;
        case 'D':
          return Type.Delete;
        case 'S':
          return Type.SoftClip;
        case 'H':
          return Type.HardClip;
        case 'N':
          return Type.Skip;
        case 'P':
          return Type.Pad;
        default:
          throw new IllegalArgumentException("Unrecognized alignment operation symbol " + sym);
      }
    }
  }

  protected static final Pattern CigarElementPattern = Pattern.compile("(\\d+)([MIDNSHP])");

  private Type op;
  private int len;

  public AlignOp(Type op, int len) {
    this.op = op;
    this.len = len;
  }

  public AlignOp(char op, int len) {
    if (len <= 0)
      throw new IllegalArgumentException("Op length must be > 0");
    this.op = Type.fromSymbol(op);
    this.len = len;
  }

  public AlignOp.Type getType() { return op; }
  public int getLen() { return len; }

  public boolean equals(Object other)
  {
    if (other instanceof AlignOp)
    {
      AlignOp otherAlignment = (AlignOp) other;
      return otherAlignment.op == this.op && otherAlignment.len == this.len;
    }
    else
      return false;
  }

  public String toString() { return "(" + op + "," + len + ")"; }

  /**
   * Compute the SAM-style CIGAR string for the given alignment.
   */
  public static String cigarStr(List<AlignOp> alignment)
  {
    if (alignment == null || alignment.isEmpty())
      return "*";

    StringBuilder builder = new StringBuilder(50);

    for (AlignOp op: alignment)
      builder.append(op.getLen()).append(op.getType().getSymbol());

    return builder.toString();
  }

  public static List<AlignOp> scanCigar(String cigar)
  {
    if ("*".equals(cigar))
      return new ArrayList<AlignOp>(0);
    else
    {
      ArrayList<AlignOp> result = new ArrayList<AlignOp>(5);
      Matcher m = CigarElementPattern.matcher(cigar);

      int lastPositionMatched = 0;
      while (m.find())
      {
        result.add( new AlignOp(AlignOp.Type.fromSymbol(m.group(2)), Integer.parseInt(m.group(1))) );
        lastPositionMatched = m.end();
      }

      if (lastPositionMatched < cigar.length())
        throw new IllegalArgumentException("Invalid CIGAR pattern " + cigar);

      if (result.isEmpty())
        throw new IllegalArgumentException("Unable to parse any alignments from CIGAR pattern " + cigar);

      return result;
    }
  }
}
