import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class SuffixTree {

	//Maps the path label of a leaf to a leaf
	private Map<String, Node> nodeMap;
	//All of the strings that have been added to the generalized suffix tree
	private List<String> stringsAdded;
	private Node root;
	//Current string we are building the suffix tree out of
	private String string;
	//Used to ensure internal Node's labels are unique
	private int internalNodeCounter;
	
	public static void main(String[] args)
	{
		SuffixTree sTree = new SuffixTree("ATAT");
		sTree.addString("TATT");
		Node root = sTree.getRoot();
		System.out.println(lce(sTree, "ATAT", 0, "TATT", 1));
	}
	
	public SuffixTree(String string)
	{	
		internalNodeCounter = 0;
		this.string = string;
		nodeMap = new HashMap<String, Node>();
		stringsAdded = new ArrayList<String>();
		//$ at end to ensure the last letter is not found in string
		constructTree(string);	
	}
	
	public void constructTree(String string)
	{
		root = new Node(new Pair<String, Integer>("r", -1), "");
		root.setFather(null);
		addString(string);	
	}
	
	public void addString(String string)
	{
		stringsAdded.add(string);
		this.string = string.concat("$");
		for (int i = 0; i < this.string.length(); i++)
		{
			addSuffix(i, root, this.string.subSequence(i, this.string.length()));
		}
	}
	
	public static String lce(SuffixTree tree, String s, int sPos, String t, int tPos)
	{
		Node sNode = tree.getLeafFromString(s.substring(sPos).concat("$"));
		Node tNode = tree.getLeafFromString(t.substring(tPos).concat("$"));
		
		if (sNode != null && tNode != null)
		{
			return lca(tree, sNode, tNode).getPathLabel();
		}
		return null;
	}
	
	public Leaf getLeafFromString(String path)
	{
		if (path != null)
		{
			return (Leaf)nodeMap.get(path);
		}
		return null;
	}
	
	public static Node lca(SuffixTree tree, Node first, Node second)
	{
		List<Node> firstAncestors = first.getAncestors();
		List<Node> secondAncestors = second.getAncestors();
		
		for (int i = 0; i < firstAncestors.size(); i++)
		{
			if (secondAncestors.contains(firstAncestors.get(i)))
			{
				return firstAncestors.get(i);
			}
		}
		return null;				
	}
	
	public void addSuffix(int startSuffix, Node father, CharSequence suffix)
	{
		if (suffix.length() == 0)
		{
			return;
		}
		List<Node> nodes = father.getChildren();
		//CharSequence suffix = string.substring(startSuffix);
		
		for (int i = 0; i < nodes.size(); i++)
		{
			Node child = nodes.get(i);
			Pair<Integer, Node> match = child.findPath(suffix);
			//First character of suffix was found
			if (match != null)
			{
				//Matched the entire edge -> Must do same thing on child
				if (match.first == child.getEdge().length())
				{
					addSuffix(startSuffix, child, suffix.subSequence(match.first, suffix.length()));
				}
				//Entire edge was not matched -> Must split and make internal Node
				else
				{
					splitEdge(father, suffix, startSuffix, match);
				}
				return;
			}
		}
		char left;
		if (startSuffix == 0)
		{
			left = '+';
		}
		else
		{
			left = string.charAt(startSuffix - 1);
		}
		//First character of suffix not found -> Add a new leaf to the node
		Node leaf = new Leaf(new Pair<String, Integer>(string, startSuffix), suffix, left);
		leaf.setFather(father);
		father.addChild(leaf);
		nodeMap.put(leaf.getPathLabel(), leaf);		
	}
	
	public void splitEdge(Node father, CharSequence suffix, int position, Pair<Integer, Node> match)
	{
		char left;
		if (position == 0)
		{
			left = '+';
		}
		else
		{
			left = string.charAt(position - 1);
		}
		father.removeChild(match.second);					
		//Create the new leaf
		Node leaf = new Leaf(
				new Pair<String, Integer>(string, position),
				suffix.subSequence(match.first, suffix.length()),
				left);					
		//Make the new internal node
		Node internal = new Node(
				new Pair<String, Integer>("Internal", internalNodeCounter++), 
				suffix.subSequence(0, match.first));	
		//Update old edge to coincide with the new internal node
		match.second.setEdge(
				match.second.getEdge().subSequence(match.first, match.second.getEdge().length()));

		internal.setFather(father);
		match.second.setFather(internal);
		leaf.setFather(internal);
		
		internal.addChild(match.second);
		internal.addChild(leaf);
		father.addChild(internal);
		
		/*if (match.second instanceof Leaf)
		{
			String path = match.second.getPathLabel();
			nodeMap.put(path, match.second);
		}*/
		
		nodeMap.put(leaf.getPathLabel(), leaf);
	}
	
	public Node getRoot()
	{
		return root;
	}
	
}
