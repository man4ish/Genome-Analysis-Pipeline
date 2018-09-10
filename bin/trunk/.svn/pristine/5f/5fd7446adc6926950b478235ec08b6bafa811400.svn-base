
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.math.BigInteger;



/**
 * This program displays factorials as the user enters values interactively
 */
public class FactQuoter {
  public static void main(String[] args) throws IOException {
    // This is how we set things up to read lines of text from the user.
    BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
    // Loop forever
    for (;;) {
      // Display a prompt to the user
      System.out.print("FactQuoter> ");
      // Read a line from the user
      String line = in.readLine();
      // If we reach the end-of-file,
      // or if the user types "quit", then quit
      if ((line == null) || line.equals("quit"))
        break;
      // Try to parse the line, and compute and print the factorial
      try {
        int x = Integer.parseInt(line);
        System.out.println(x + "! = " + Factorial4.factorial(x));
      }
      // If anything goes wrong, display a generic error message
      catch (Exception e) {
        System.out.println("Invalid Input");
      }
    }
  }
}
/**
 * This version of the program uses arbitrary precision integers, so it does not
 * have an upper-bound on the values it can compute. It uses an ArrayList object
 * to cache computed values instead of a fixed-size array. An ArrayList is like
 * an array, but can grow to any size. The factorial() method is declared
 * "synchronized" so that it can be safely used in multi-threaded programs. Look
 * up java.math.BigInteger and java.util.ArrayList while studying this class.
 * Prior to Java 1.2, use Vector instead of ArrayList
 */

class Factorial4 {
  protected static ArrayList table = new ArrayList(); // create cache
  static { // Initialize the first element of the cache with !0 = 1.
    table.add(BigInteger.valueOf(1));
  }

  /** The factorial() method, using BigIntegers cached in a ArrayList */
  public static synchronized BigInteger factorial(int x) {
    if (x < 0)
      throw new IllegalArgumentException("x must be non-negative.");
    for (int size = table.size(); size <= x; size++) {
      BigInteger lastfact = (BigInteger) table.get(size - 1);
      BigInteger nextfact = lastfact.multiply(BigInteger.valueOf(size));
      table.add(nextfact);
    }
    return (BigInteger) table.get(x);
  }

  /**
   * A simple main() method that we can use as a standalone test program for
   * our factorial() method.
   */
  /*public static void main(String[] args) {
    for (int i = 0; i <= 50; i++)
      System.out.println(i + "! = " + factorial(i));
  }
  */
}

