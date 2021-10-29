#include "FamilyIdentifierHandler.h"
#include "SmilesHandler.h"
#include <cmath>

namespace combi_ff{

namespace cnv{

void FamilyIdentifierHandler::Run() {
	cnv::Handler::Run();

	//determine the different smiles strings from the input_list
	for(auto fmi_it = input_list.begin(); fmi_it != input_list.end(); ++fmi_it) {
		for(size_t i = 0; i < fmi_it->second.size(); i++) {
			if((fmi_it->second)[i] == ' ' || (fmi_it->second)[i] == ',') {
				auto next = ++fmi_it;
				fmi_it--;
				input_list.insert(next, {fmi_it->first, fmi_it->second.substr(i + 1, fmi_it->second.size() - i)});
				fmi_it->second = fmi_it->second.substr(0, i);
			}
		}
	}

	for(auto&& fmi : input_list) {
		if(fmi.second.size())
			FamilyIdentifierHandler::ConvertFamilyIdentifier(fmi.second);
	}
}

void FamilyIdentifierHandler::PrintFirstLine() {
	cnv::Handler::PrintFirstLine("# familyIdentifier ");
}

void FamilyIdentifierHandler::ConvertFamilyIdentifier(const std::string& fmi) {
	std::string family_code = fmi.substr(0, 4);
	bool found(false);

	for(auto && family : fie_file_names) {
		if(family.find(family_code) != std::string::npos) {
			found = true;
			FamilyIdentifierHandler::SearchFileForIsomer(fmi, family);
			break;
		}
	}

	if(!found) {
		std::cout << "!Error: no family with the code " << family_code << " for compound " << fmi << " found in any of the family isomer enumeration files ";

		for(auto && ffn : fie_file_names)
			std::cout << ffn << " ";

		std::cout << '\n';
		exit(-1);
	}
}

void FamilyIdentifierHandler::SearchFileForIsomer(const std::string& fmi, const std::string& family_file_name) {
	std::ifstream family_file(family_file_name);

	if(!family_file.is_open()) {
		std::cout << "!Error: couldn't open family file " << family_file_name << " for compound " << fmi << "\n";
		exit(-1);
	}

	bool found(false);
	std::string line("");

	while(std::getline(family_file, line)) {
		if(line[0] != '#' && line.find(fmi) != std::string::npos) {
			found = true;
			break;
		}
	}

	if(!found) {
		std::cout << "!Error: couldn't find compound " << fmi << " in " << family_file_name << "\n";
		exit(-1);
	}

	else {
		family_file.seekg((int)family_file.tellg() - (line.size() + 1));
		std::string code, formula, smiles_orig;
		family_file >> code >> formula >> smiles_orig;
		//assert(code == fmi);
		cnv::AdjacencyMatrix A;
		std::string smiles_canon;

		if(print_options[cnv::print_mass] || print_options[cnv::print_stack] || print_options[cnv::print_matrix] || print_options[cnv::print_n_unsaturations] || print_options[cnv::print_n_multiple_bonds] ||
		        print_options[cnv::print_n_double_bonds] || print_options[cnv::print_n_aromatic_bonds] || print_options[cnv::print_n_triple_bonds] || print_options[cnv::print_n_cycles] || print_options[cnv::print_n_quadruple_bonds]){
			combi_ff::cnv::SmilesHandler(output_file, fie_file_names, column_width, print_options, input_list).ConvertSmiles(smiles_orig, smiles_canon, A);
			//assert(smiles_orig == smiles_canon);
		}

		FamilyIdentifierHandler::PrintOutput(fmi, formula, smiles_orig, A);
	}

	family_file.close();
}



void FamilyIdentifierHandler::PrintOutput(const std::string& fmi, const std::string& formula, const std::string& smiles, AdjacencyMatrix& A) {
	size_t nUnsat(0), nMB(0), nDB(0), nAB(0), nTB(0), nSB(0), nQB(0), nCyc(0);

	if(print_options[cnv::print_n_double_bonds] || print_options[cnv::print_n_aromatic_bonds] || print_options[cnv::print_n_triple_bonds] || print_options[cnv::print_n_unsaturations] || print_options[cnv::print_n_cycles] ||
	        print_options[cnv::print_n_multiple_bonds]) {
		A.GetNumMultipleBonds(nSB, nDB, nTB, nQB, nAB);
		double DoU_(0);
		//assert(!(nAB % 2));
		nMB = nDB + nTB + nAB / 2;

		for(auto && atom : A.GetAtomVector())
			DoU_ += (double)atom.GetDegree() - 2.;

		nUnsat = size_t(floor((DoU_ / 2.) + 1));
		nCyc = nUnsat - nMB;
	}

	std::ostream* out;

	if(output_file.is_open())
		out = &output_file;

	else
		out = &std::cout;

	*out << std::setw(column_width) << std::left << fmi << " ";

	if(print_options[cnv::print_canon_smiles])
		*out << std::setw(column_width) << std::left << smiles << " ";

	if(print_options[cnv::print_canon_formula])
		*out << std::setw(column_width) << std::left << formula << " ";

	if(print_options[cnv::print_mass]) {
		double mass = 0;

		for(auto && a : A.GetAtomVector())
			mass += a.GetMass();

		*out << std::setw(column_width) << std::left << mass << " ";
	}

	if(print_options[cnv::print_n_unsaturations])
		*out << std::setw(column_width) << std::left << nUnsat << " ";

	if(print_options[cnv::print_n_multiple_bonds])
		*out << std::setw(column_width) << std::left << nMB << " ";

	if(print_options[cnv::print_n_double_bonds])
		*out << std::setw(column_width) << std::left << nDB << " ";

	if(print_options[cnv::print_n_aromatic_bonds])
		*out << std::setw(column_width) << std::left << nAB << " ";

	if(print_options[cnv::print_n_triple_bonds])
		*out << std::setw(column_width) << std::left << nTB << " ";

	if(print_options[cnv::print_n_quadruple_bonds])
		*out << std::setw(column_width) << std::left << nQB << " ";

	if(print_options[cnv::print_n_cycles])
		*out << std::setw(column_width) << std::left << nCyc << " ";

	if(print_options[cnv::print_canon_atom_vector])
		*out << A.GetAtomVector() << " ";

	if(print_options[cnv::print_stack])
		*out << A.GetStack();

	if(print_options[cnv::print_matrix]) {
		*out << '\n';
		A.PrintToFile(*out);
	}

	*out << '\n';
}


} //namespace cnv

} //namespace combi_ff