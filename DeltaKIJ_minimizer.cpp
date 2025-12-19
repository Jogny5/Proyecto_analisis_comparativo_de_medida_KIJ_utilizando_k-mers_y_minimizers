#include "HLL_Minimizers.cpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <sstream>
using namespace std;

ofstream verbose;
ofstream verbose_min;

unordered_map<string, double> delta_cache;       // normales
unordered_map<string, double> delta_cache_min;   // minimizers

uint W = 10;                                     // ventana de minimizers
const size_t BLOCK_SIZE = 8;                     // tamaño del bloque 2, 4, 8
using HLLVec = vector<HyperLogLog>; 

HLLVec make_hll_vec() {
    HLLVec v;
    for (uint k = 3; k <= 24; k++)
        v.emplace_back(14);
    return v;
}

// Calcula el Hyperloglog para los archivos
HLLVec build_hlls_for_file(const string& file) {
    HLLVec hlls = make_hll_vec();

    for (uint k = 3; k <= 24; k++) {
        process_kmers_hll(hlls[k-2], file, k);
    }
    return hlls;
}

// Calcula el Hyperloglog para los archivos usando minimizers
HLLVec build_hlls_for_file_min(const string& file) {
    HLLVec hlls = make_hll_vec(); 

    for (uint k = 3; k <= 24; k++) {
        process_kmers_minimizers(hlls[k-3], file, k, W);
    }
    return hlls;
}

// Calcula el delta con Hyperloglog
double compute_delta_from_hll(HyperLogLog &hll) {
    double max_delta = 0.0;

    for (uint k = 3; k <= 24; k++) {
        double est = hll.estimate();
        double d = est / (double)k;
        max_delta = max(max_delta, d);
    }

    return max_delta;
}
// Calcula el delta con Hyperloglog usando minimizers
double compute_delta_from_hlls_min(const HLLVec& hlls) {
    double max_delta = 0.0;

    for (uint k = 3; k <= 24; k++) {
        double est = hlls[k-3].estimate();
        double d = est / (double)k;
        max_delta = max(max_delta, d);
    }

    return max_delta;
}

// función que usa merge para unir archivos
void merge_hlls(HLLVec& base, const HLLVec& other) {
    for (size_t i = 0; i < base.size(); i++) {
        base[i].merge_from(other[i]);
    }
}




void dual_print(const string &s) {
    cout << s;
    verbose << s;
}

void dual_print_min(const string &s) {
    cout << s;
    verbose_min << s;
}

void load_delta_cache() {
    ifstream in("delta_cache.txt");
    string file; double val;
    if (in.is_open()) {
        while (in >> file >> val) delta_cache[file] = val;
        in.close();
    }

    ifstream in2("delta_cache_min.txt");
    if (in2.is_open()) {
        while (in2 >> file >> val) delta_cache_min[file] = val;
        in2.close();
    }
}

// Funcion para evitar volver a calcular deltas individuales para cada ejecución
void save_delta_cache() {
    ofstream out("delta_cache.txt");
    for (auto &p : delta_cache) out << p.first << " " << p.second << "\n";
    out.close();

    ofstream out2("delta_cache_min.txt");
    for (auto &p : delta_cache_min) out2 << p.first << " " << p.second << "\n";
    out2.close();
}

// Calculo de delta para archivo sin minimizers
double compute_delta_normal_raw(const string& filename) {
    dual_print("\n=== ARCHIVO: " + filename + " ===\n");

    double max_delta = 0;

    for (uint k = 2; k <= 31; k++) {
        HyperLogLog hll(14);
        process_kmers_hll(hll, filename, k);

        double est = hll.estimate();
        double d = est / (double)k;

        std::ostringstream ss;
        ss << "  k=" << k
           << "  H(k)=" << (uint64_t)est
           << "   H(k)/k=" << d << "\n";
        dual_print(ss.str());

        max_delta = max(max_delta, d);
    }

    return max_delta;
}


// Calculo de delta para archivo con minimizers
double compute_delta_minimizers_raw(const string& filename) {
    dual_print_min("\n=== ARCHIVO (MIN): " + filename + " ===\n");

    double max_delta = 0;

    for (uint k = 2; k <= 31; k++) {
        HyperLogLog hll(14);
        process_kmers_minimizers(hll, filename, k, W);

        double est = hll.estimate();
        double d = est / (double)k;

        std::ostringstream ss;
        ss << "[MIN] k=" << k
           << "  H(k)=" << (uint64_t)est
           << "   H(k)/k=" << d << "\n";
        dual_print_min(ss.str());

        max_delta = max(max_delta, d);
    }

    return max_delta;
}

// Busca si ya está calculado el delta de un archivo
double get_delta_normal(const string& filename) {
    if (delta_cache.count(filename)) {
        dual_print("[CACHE] delta(" + filename + ") = "
                   + to_string(delta_cache[filename]) + "\n");
        return delta_cache[filename];
    }

    double d = compute_delta_normal_raw(filename);
    delta_cache[filename] = d;
    save_delta_cache();
    return d;
}

double get_delta_minimizers(const string& filename) {
    if (delta_cache_min.count(filename)) {
        dual_print_min("[CACHE] delta_min(" + filename + ") = "
                       + to_string(delta_cache_min[filename]) + "\n");
        return delta_cache_min[filename];
    }

    double d = compute_delta_minimizers_raw(filename);
    delta_cache_min[filename] = d;
    save_delta_cache();
    return d;
}

// Función para la unión de los archivos
double compute_delta_union(const string& A, const string& B) {

    dual_print("\n=== UNION: " + A + " ∪ " + B + " ===\n");

    double max_delta = 0;

    for (uint k = 3; k <= 24; k++) {
        HyperLogLog hA(14), hB(14);
        process_kmers_hll(hA, A, k);
        process_kmers_hll(hB, B, k);
        hA.merge_from(hB);

        double est = hA.estimate();
        double d = est / (double)k;

        std::ostringstream ss;
        ss << "  k=" << k << "  H(k)=" << (uint64_t)est
           << "   H(k)/k=" << d << "\n";
        dual_print(ss.str());

        max_delta = max(max_delta, d);
    }

    return max_delta;
}

// Función para la unión de los archivos con minimizers
double compute_delta_union_min(const string& A, const string& B) {

    dual_print_min("\n=== UNION (MIN): " + A + " ∪ " + B + " ===\n");

    double max_delta = 0;

    for (uint k = 2; k <= 31; k++) {
        HyperLogLog hA(14), hB(14);
        process_kmers_minimizers(hA, A, k, W);
        process_kmers_minimizers(hB, B, k, W);
        hA.merge_from(hB);

        double est = hA.estimate();
        double d = est / (double)k;

        std::ostringstream ss;
        ss << "[MIN] k=" << k << "  H(k)=" << (uint64_t)est
           << "   H(k)/k=" << d << "\n";
        dual_print_min(ss.str());

        max_delta = max(max_delta, d);
    }

    return max_delta;
}

// Calcula el kij 
double compute_kij(double dA, double dB, double dU) {
    double kij = (dA + dB - dU) / dU;
    if (kij < 0) kij = 0;
    if (kij > 1) kij = 1;
    return kij;
}

// Imprime los resultados y los guarda en los archivos
void print_and_save_kij(
    const string& label,
    double dA,
    double dB,
    double dU,
    double kij,
    ofstream& out,
    ofstream& verbose
) {
    std::ostringstream ss;

    ss << "\n=== " << label << " ===\n";
    ss << "delta(A): " << dA << "\n";
    ss << "delta(B): " << dB << "\n";
    ss << "delta(A ∪ B): " << dU << "\n";
    ss << "Kij: " << kij << "\n\n";

    cout << ss.str();
    verbose << ss.str();

    out << "delta(A) = " << (uint64_t)dA << "\n";
    out << "delta(B) = " << (uint64_t)dB << "\n";
    out << "delta(A ∪ B) = " << (uint64_t)dU << "\n";
    out << "Kij = " << kij << "\n\n";

    out.flush();
    verbose.flush();
}



int main(int argc, char* argv[]) {

    if (argc < 2) {
        cerr << "Uso: " << argv[0] << " lista_archivos.txt\n";
        return 1;
    }

    load_delta_cache();

    string list = argv[1];
    ifstream in(list);

    vector<string> files;
    string line;
    while (getline(in, line)) {
        line.erase(line.find_last_not_of(" \r\n\t") + 1);
        if (!line.empty()) files.push_back(line);
    }
    in.close();

    //archivos para guardar los resultados
    ofstream out("kij_results.txt");
    ofstream out_min("kij_results_minimizers.txt");

    verbose.open("kij_verbose_output.txt");
    verbose_min.open("kij_verbose_output_minimizers.txt");

    // Union de bloques sin minimizers

    // for (size_t b = 0; b < files.size(); b += BLOCK_SIZE) {

    //     size_t end = min(b + BLOCK_SIZE, files.size());
    //     if (end - b < 2) break;

    //     dual_print("\n========================================\n");
    //     dual_print(" NUEVO BLOQUE\n");
    //     dual_print("========================================\n");

    //     // HLL acumulados por k
    //     vector<HyperLogLog> hll_union(32, HyperLogLog(14));

    //     // ---------- Inicializar con A ----------
    //     string A = files[b];
    //     for (uint k = 2; k <= 31; k++)
    //         process_kmers_hll(hll_union[k], A, k);

    //     double dA = get_delta_normal(A);

    //     // ---------- B ----------
    //     string B = files[b + 1];
    //     for (uint k = 2; k <= 31; k++)
    //         process_kmers_hll(hll_union[k], B, k);

    //     double dB = get_delta_normal(B);

    //     double dAB = 0.0;
    //     for (uint k = 2; k <= 31; k++)
    //         dAB = max(dAB, hll_union[k].estimate() / (double)k);

    //     double kij = compute_kij(dA, dB, dAB);

    //     print_and_save_kij(
    //         "PAR: " + A + " y " + B,
    //         dA, dB, dAB, kij,
    //         out, verbose
    //     );

    //     double dPrev = dAB;
    //     string label = "(" + A + " ∪ " + B + ")";

    //     // ---------- Incremental ----------
    //     for (size_t i = b + 2; i < end; i++) {

    //         string C = files[i];
    //         double dC = get_delta_normal(C);

    //         for (uint k = 2; k <= 31; k++)
    //             process_kmers_hll(hll_union[k], C, k);

    //         double dNew = 0.0;
    //         for (uint k = 2; k <= 31; k++)
    //             dNew = max(dNew, hll_union[k].estimate() / (double)k);

    //         double kij_inc = compute_kij(dPrev, dC, dNew);

    //         print_and_save_kij(
    //             "PAR: " + label + " y " + C,
    //             dPrev, dC, dNew, kij_inc,
    //             out, verbose
    //         );

    //         dPrev = dNew;
    //         label = "(" + label + " ∪ " + C + ")";
    //     }
    // }

    // -------- Para las ejecuciones con minimizers--------

    for (size_t b = 0; b < files.size(); b += BLOCK_SIZE) {

        size_t end = min(b + BLOCK_SIZE, files.size());
        if (end - b < 2) break;

        dual_print_min("\n========================================\n");
        dual_print_min(" NUEVO BLOQUE (MINIMIZERS)\n");
        dual_print_min("========================================\n");

        string A = files[b];
        HLLVec hll_union_min = build_hlls_for_file_min(A);
        double dAmin = get_delta_minimizers(A);

        string B = files[b + 1];
        HLLVec hll_B_min = build_hlls_for_file_min(B);
        double dBmin = get_delta_minimizers(B);

        merge_hlls(hll_union_min, hll_B_min);

        double dABmin = compute_delta_from_hlls_min(hll_union_min);
        double kijAB = compute_kij(dAmin, dBmin, dABmin);

        print_and_save_kij(
            "PAR (MIN): " + A + " y " + B,
            dAmin, dBmin, dABmin, kijAB,
            out_min, verbose_min
        );

        double dPrev = dABmin;
        string label = "(" + A + " ∪ " + B + ")";

        for (size_t i = b + 2; i < end; i++) {

            string C = files[i];
            double dCmin = get_delta_minimizers(C);

            HLLVec hll_C_min = build_hlls_for_file_min(C);
            merge_hlls(hll_union_min, hll_C_min);

            double dNew = compute_delta_from_hlls_min(hll_union_min);
            double kij_inc = compute_kij(dPrev, dCmin, dNew);

            print_and_save_kij(
                "PAR (MIN): " + label + " y " + C,
                dPrev, dCmin, dNew, kij_inc,
                out_min, verbose_min
            );

            dPrev = dNew;
            label = "(" + label + " ∪ " + C + ")";
        }
    }

    verbose.close();
    verbose_min.close();
    out.close();
    out_min.close();

    cout << "\nListo. Resultados:\n";
    cout << " - kij_results.txt\n";
    cout << " - kij_verbose_output.txt\n";
    cout << " - kij_results_minimizers.txt\n";
    cout << " - kij_verbose_output_minimizers.txt\n";
    return 0;
}
