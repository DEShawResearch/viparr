#include "system_to_dot.hxx"

using namespace desres::msys;

namespace desres{ namespace viparr {

    void SystemToDot(SystemPtr sys, std::ostream& dotfile, Id residue_id) {
        dotfile << "graph topology {\n";
        dotfile << "\tnode [fontsize=11];\n";
        IdList atoms = sys->atoms();
        for (unsigned i = 0; i < atoms.size(); ++i) {
            Id atom = atoms[i];
            int anum = sys->atom(atom).atomic_number;
            std::string name = sys->atom(atom).name;
            if (!bad(residue_id) && sys->atom(atom).residue != residue_id) {
                anum = -1;
            }
            switch (anum) {
                case -1: // External
                    dotfile << "\t\"" << atom << "\" [shape=tripleoctagon, fillcolor=yellow, style=filled, label=\"" << name << "\"];\n";
                    break;
                case 0: // Virtual
                    dotfile << "\t\"" << atom << "\" [shape=rect, width=0, height=0, margin=0, label=\"" << name << "\\nvsite\"];\n";
                    break;
                case 1: // H
                    dotfile << "\t\"" << atom << "\" [shape=circle, width=0, height=0, margin=0, label=\"" << name << "\\nH(" << atom << ")\"];\n";
                    break;
                case 6: // C
                    dotfile << "\t\"" << atom << "\" [shape=circle, fillcolor=beige, style=filled, label=\"" << name << "\\nC(" << atom << ")\"];\n";
                    break;
                case 7: // N
                    dotfile << "\t\"" << atom << "\" [shape=house, fillcolor=coral, style=filled, label=\"" << name << "\\nN(" << atom << ")\"];\n";
                    break;
                case 8: // O
                    dotfile << "\t\"" << atom << "\" [shape=circle, fillcolor=cyan, style=filled, label=\"" << name << "\\nO(" << atom << ")\"];\n";
                    break;
                case 15: // P
                    dotfile << "\t\"" << atom << "\" [shape=octagon, fillcolor=red, style=filled, label=\"" << name << "\\nP(" << atom << ")\"];\n";
                    break;
                case 16: // S
                    dotfile << "\t\"" << atom << "\" [shape=hexagon, fillcolor=gold, style=filled, label=\"" << name << "\\nS(" << atom << ")\"];\n";
                    break;
                default: // Other
                    dotfile << "\t\"" << atom << "\" [shape=tripleoctagon, fillcolor=yellow, style=filled, label=\"" << name << "\\n" << anum << "(" << atom << ")\"];\n";
                    break;
            }
        }
        dotfile << "\n";
        IdList bonds = sys->bonds();
        for (unsigned i = 0; i < bonds.size(); ++i)
            dotfile << "\t\"" << sys->bond(bonds[i]).i << "\" -- \"" << sys->bond(bonds[i]).j << "\";\n";
        dotfile << "}" << std::endl;
    }

}}
