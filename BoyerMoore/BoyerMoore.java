import java.nio.CharBuffer;
import java.util.HashMap;
import java.util.Map;



public class BoyerMoore {
	
	private CharSequence sub = null;
	//private int substringLength;
	
	private CharSequence alpha = null;
	//private int alphaLength;
	
	private Map<Character, Integer> badCharRule;
	
	private int[] prefixValues;
	
	public BoyerMoore(CharSequence alphabet, CharSequence substring)
	{
		alpha = alphabet;
		sub = substring;
		/*substringLength = sub.length();
		alphaLength = alpha.length();*/
		
		prefixValues = prefixFunction(sub);
		badCharRule = computeBadCharRule(alpha, sub);
		computeGoodSuffixShift(sub);
	}
	
	private int BM(CharSequence str, CharSequence substring, CharSequence alphabet)
	{
		Map<Character, Integer> R = computeBadCharRule(substring, alphabet);
		int[] GS = computeGoodSuffixShift(substring);
		int n = str.length();
		int m = substring.length();
		int s = 0;
		
		while (s <= n - m)
		{
			int j = m;
			while (j > 0 && substring.charAt(j) == str.charAt(s + j))
			{
				j--;
			}
			if (j == 0)
			{
				return s;
				//s = s + GS[0];
			}
			else
			{
				s = Math.max(GS[j], j - R.get(str.charAt(s+j)));
			}
		}
		return -1;
	}
	
	private Map<Character, Integer> computeBadCharRule(CharSequence substring, CharSequence alphabet)
	{
		if (!(alphabet != null && substring != null))
		{
			return null;
		}
		Map<Character, Integer> badRule = new HashMap<Character, Integer>(); 
		
		for (int i = 0; i < alphabet.length(); i++)
		{
			badRule.put(alphabet.charAt(i), 0);
		}
		for (int i = 0; i < substring.length(); i++)
		{
			badRule.put(substring.charAt(i), i);
		}
		return badRule;
	}
	
	private int[] computeGoodSuffixShift(CharSequence substring)
	{
		if (substring == null)
		{
			return null;
		}
		int length = substring.length();
		int[] goodRule = new int[length];
		
		int[] prefixForward = prefixFunction(substring);
		String reversedSubstring = new StringBuilder(substring).reverse().toString();
		int[] prefixReverse = prefixFunction(reversedSubstring);
		
		for (int i = 0; i < length; i++)
		{
			goodRule[i] = length-1 - prefixForward[length-1];
		}
		for(int i = 0; i < length; i++)
		{
			int j = length - 1 - prefixReverse[i];
			if (goodRule[j] > (i - prefixReverse[i]))
			{
				goodRule[j] = i - prefixReverse[i];
			}
		}
		return goodRule;
	}
	
	private int[] prefixFunction(CharSequence substring)
	{
		if (substring == null)
		{
			return null;
		}
		
		int[] piValues = new int[substring.length()];
		piValues[0] = 0;
		int k = 0;
		for (int i = 1; i < substring.length(); i++)
		{
			while (k > 0 && substring.charAt(k) != substring.charAt(i))
			{
				k = piValues[k-1];
			}
			if (substring.charAt(k) == substring.charAt(i))
			{
				k++;
			}
			piValues[i] = k;
		}
		return piValues;
	}
	
	public static void main(String[] args)
	{
		BoyerMoore bm = new BoyerMoore("GCAT", "GTACTTACTTAC");
	}

}
