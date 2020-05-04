
#ifndef TREE_H
#define TREE_H

#include <string>

#include "Params.h++"


/**********************************/
/******      CLASS TREE      ******/
/**********************************/

class Tree
{
	public :
	// constructors
		Tree(par::Params const &params);
		Tree();

	// friend functions and overload
//		Tree& operator=(const Tree& tree);
		void printTree(std::ostream &os) const;

	// accessors
		// std::string getName() const;
		//
		// double getOverG() const;
		// double getUnderG() const;
		// double getOverMu() const;
		// double getUnderMu() const;
		//
		// double getFecundity() const;
		// double getPhi() const;
		// double getMaxDiameter() const;
		//
		// double getAlpha20() const;
		// double getBeta() const;
		// double getAlpha() const;

	// others
		bool compareName(Tree* tree_1, Tree* tree_2) const;

	private :
		std::string m_treeName;
		double m_overG, m_underG, m_overMu, m_underMu, m_phi, m_fecundity, m_maxDiameter, m_alpha20, m_beta, m_alpha;
		// WATCH OUT, alpha != alpha20, cf Purves 2007/2008.
};

std::ostream& operator<<(std::ostream& os, const Tree& tree);

#endif
