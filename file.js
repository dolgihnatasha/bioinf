let fs = require('fs');

let useful = {
    "codonToAminoRNA": {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M", "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S", "CCU": "P", "CCC": "P",
        "CCA": "P", "CCG": "P", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCU": "A", "GCC": "A",
        "GCA": "A", "GCG": "A", "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*", "UGA": "*", "CAU": "H",
        "CAC": "H", "CAA": "Q", "CAG": "Q", "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAU": "D",
        "GAC": "D", "GAA": "E", "GAG": "E", "UGU": "C", "UGC": "C", "UGG": "W", "CGU": "R", "CGC": "R",
        "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R", "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    },
    "codonToAminoDNA": {
        'TTT': 'F',     'CTT': 'L',     'ATT': 'I',     'GTT': 'V',
        'TTC': 'F',     'CTC': 'L',     'ATC': 'I',     'GTC': 'V',
        'TTA': 'L',     'CTA': 'L',     'ATA': 'I',     'GTA': 'V',
        'TTG': 'L',     'CTG': 'L',     'ATG': 'M',     'GTG': 'V',
        'TCT': 'S',     'CCT': 'P',     'ACT': 'T',     'GCT': 'A',
        'TCC': 'S',     'CCC': 'P',     'ACC': 'T',     'GCC': 'A',
        'TCA': 'S',     'CCA': 'P',     'ACA': 'T',     'GCA': 'A',
        'TCG': 'S',     'CCG': 'P',     'ACG': 'T',     'GCG': 'A',
        'TAT': 'Y',     'CAT': 'H',     'AAT': 'N',     'GAT': 'D',
        'TAC': 'Y',     'CAC': 'H',     'AAC': 'N',     'GAC': 'D',
        'TAA': '*',     'CAA': 'Q',     'AAA': 'K',     'GAA': 'E',
        'TAG': '*',     'CAG': 'Q',     'AAG': 'K',     'GAG': 'E',
        'TGT': 'C',     'CGT': 'R',     'AGT': 'S',     'GGT': 'G',
        'TGC': 'C',     'CGC': 'R',     'AGC': 'S',     'GGC': 'G',
        'TGA': '*',     'CGA': 'R',     'AGA': 'R',     'GGA': 'G',
        'TGG': 'W',     'CGG': 'R',     'AGG': 'R',     'GGG': 'G'
    },
    "aminoToCodon": {
        "F": "UUC", "L": "CUG", "I": "AUC", "M": "AUG", "V": "GUG", "S": "AGC",
        "P": "CCC", "T": "ACC", "A": "GCC", "Y": "AUC", "*": "UGA", "H": "CAC",
        "Q": "CAG", "N": "AAC", "K": "AAG", "D": "GAC", "E": "GAG", "C": "UGC",
        "R": "AGA", "G": "GGC"
    },
    "DNAToRNA": {
        "A": "U",
        "T": "A",
        "C": "G",
        "G": "C"
    },
    "DNAToDNA": {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    },
    "RNAToDNA": {
        "U": "A",
        "A": "T",
        "G": "C",
        "C": "G"
    },
    "RNAToRNA": {
        "U": "A",
        "A": "U",
        "G": "C",
        "C": "G"
    },
    "blosum62": {
        "-": {
            "-": 1,
            "A": -4,
            "C": -4,
            "B": -4,
            "E": -4,
            "D": -4,
            "G": -4,
            "F": -4,
            "I": -4,
            "H": -4,
            "K": -4,
            "M": -4,
            "L": -4,
            "N": -4,
            "Q": -4,
            "P": -4,
            "S": -4,
            "R": -4,
            "T": -4,
            "W": -4,
            "V": -4,
            "Y": -4,
            "X": -4,
            "Z": -4
        },
        "A": {
            "-": -4,
            "A": 4,
            "C": 0,
            "B": -2,
            "E": -1,
            "D": -2,
            "G": 0,
            "F": -2,
            "I": -1,
            "H": -2,
            "K": -1,
            "M": -1,
            "L": -1,
            "N": -2,
            "Q": -1,
            "P": -1,
            "S": 1,
            "R": -1,
            "T": 0,
            "W": -3,
            "V": 0,
            "Y": -2,
            "X": 0,
            "Z": -1
        },
        "C": {
            "-": -4,
            "A": 0,
            "C": 9,
            "B": -3,
            "E": -4,
            "D": -3,
            "G": -3,
            "F": -2,
            "I": -1,
            "H": -3,
            "K": -3,
            "M": -1,
            "L": -1,
            "N": -3,
            "Q": -3,
            "P": -3,
            "S": -1,
            "R": -3,
            "T": -1,
            "W": -2,
            "V": -1,
            "Y": -2,
            "X": -2,
            "Z": -3
        },
        "B": {
            "-": -4,
            "A": -2,
            "C": -3,
            "B": 4,
            "E": 1,
            "D": 4,
            "G": -1,
            "F": -3,
            "I": -3,
            "H": 0,
            "K": 0,
            "M": -3,
            "L": -4,
            "N": 3,
            "Q": 0,
            "P": -2,
            "S": 0,
            "R": -1,
            "T": -1,
            "W": -4,
            "V": -3,
            "Y": -3,
            "X": -1,
            "Z": 1
        },
        "E": {
            "-": -4,
            "A": -1,
            "C": -4,
            "B": 1,
            "E": 5,
            "D": 2,
            "G": -2,
            "F": -3,
            "I": -3,
            "H": 0,
            "K": 1,
            "M": -2,
            "L": -3,
            "N": 0,
            "Q": 2,
            "P": -1,
            "S": 0,
            "R": 0,
            "T": -1,
            "W": -3,
            "V": -2,
            "Y": -2,
            "X": -1,
            "Z": 4
        },
        "D": {
            "-": -4,
            "A": -2,
            "C": -3,
            "B": 4,
            "E": 2,
            "D": 6,
            "G": -1,
            "F": -3,
            "I": -3,
            "H": -1,
            "K": -1,
            "M": -3,
            "L": -4,
            "N": 1,
            "Q": 0,
            "P": -1,
            "S": 0,
            "R": -2,
            "T": -1,
            "W": -4,
            "V": -3,
            "Y": -3,
            "X": -1,
            "Z": 1
        },
        "G": {
            "-": -4,
            "A": 0,
            "C": -3,
            "B": -1,
            "E": -2,
            "D": -1,
            "G": 6,
            "F": -3,
            "I": -4,
            "H": -2,
            "K": -2,
            "M": -3,
            "L": -4,
            "N": 0,
            "Q": -2,
            "P": -2,
            "S": 0,
            "R": -2,
            "T": -2,
            "W": -2,
            "V": -3,
            "Y": -3,
            "X": -1,
            "Z": -2
        },
        "F": {
            "-": -4,
            "A": -2,
            "C": -2,
            "B": -3,
            "E": -3,
            "D": -3,
            "G": -3,
            "F": 6,
            "I": 0,
            "H": -1,
            "K": -3,
            "M": 0,
            "L": 0,
            "N": -3,
            "Q": -3,
            "P": -4,
            "S": -2,
            "R": -3,
            "T": -2,
            "W": 1,
            "V": -1,
            "Y": 3,
            "X": -1,
            "Z": -3
        },
        "I": {
            "-": -4,
            "A": -1,
            "C": -1,
            "B": -3,
            "E": -3,
            "D": -3,
            "G": -4,
            "F": 0,
            "I": 4,
            "H": -3,
            "K": -3,
            "M": 1,
            "L": 2,
            "N": -3,
            "Q": -3,
            "P": -3,
            "S": -2,
            "R": -3,
            "T": -1,
            "W": -3,
            "V": 3,
            "Y": -1,
            "X": -1,
            "Z": -3
        },
        "H": {
            "-": -4,
            "A": -2,
            "C": -3,
            "B": 0,
            "E": 0,
            "D": -1,
            "G": -2,
            "F": -1,
            "I": -3,
            "H": 8,
            "K": -1,
            "M": -2,
            "L": -3,
            "N": 1,
            "Q": 0,
            "P": -2,
            "S": -1,
            "R": 0,
            "T": -2,
            "W": -2,
            "V": -3,
            "Y": 2,
            "X": -1,
            "Z": 0
        },
        "K": {
            "-": -4,
            "A": -1,
            "C": -3,
            "B": 0,
            "E": 1,
            "D": -1,
            "G": -2,
            "F": -3,
            "I": -3,
            "H": -1,
            "K": 5,
            "M": -1,
            "L": -2,
            "N": 0,
            "Q": 1,
            "P": -1,
            "S": 0,
            "R": 2,
            "T": -1,
            "W": -3,
            "V": -2,
            "Y": -2,
            "X": -1,
            "Z": 1
        },
        "M": {
            "-": -4,
            "A": -1,
            "C": -1,
            "B": -3,
            "E": -2,
            "D": -3,
            "G": -3,
            "F": 0,
            "I": 1,
            "H": -2,
            "K": -1,
            "M": 5,
            "L": 2,
            "N": -2,
            "Q": 0,
            "P": -2,
            "S": -1,
            "R": -1,
            "T": -1,
            "W": -1,
            "V": 1,
            "Y": -1,
            "X": -1,
            "Z": -1
        },
        "L": {
            "-": -4,
            "A": -1,
            "C": -1,
            "B": -4,
            "E": -3,
            "D": -4,
            "G": -4,
            "F": 0,
            "I": 2,
            "H": -3,
            "K": -2,
            "M": 2,
            "L": 4,
            "N": -3,
            "Q": -2,
            "P": -3,
            "S": -2,
            "R": -2,
            "T": -1,
            "W": -2,
            "V": 1,
            "Y": -1,
            "X": -1,
            "Z": -3
        },
        "N": {
            "-": -4,
            "A": -2,
            "C": -3,
            "B": 3,
            "E": 0,
            "D": 1,
            "G": 0,
            "F": -3,
            "I": -3,
            "H": 1,
            "K": 0,
            "M": -2,
            "L": -3,
            "N": 6,
            "Q": 0,
            "P": -2,
            "S": 1,
            "R": 0,
            "T": 0,
            "W": -4,
            "V": -3,
            "Y": -2,
            "X": -1,
            "Z": 0
        },
        "Q": {
            "-": -4,
            "A": -1,
            "C": -3,
            "B": 0,
            "E": 2,
            "D": 0,
            "G": -2,
            "F": -3,
            "I": -3,
            "H": 0,
            "K": 1,
            "M": 0,
            "L": -2,
            "N": 0,
            "Q": 5,
            "P": -1,
            "S": 0,
            "R": 1,
            "T": -1,
            "W": -2,
            "V": -2,
            "Y": -1,
            "X": -1,
            "Z": 3
        },
        "P": {
            "-": -4,
            "A": -1,
            "C": -3,
            "B": -2,
            "E": -1,
            "D": -1,
            "G": -2,
            "F": -4,
            "I": -3,
            "H": -2,
            "K": -1,
            "M": -2,
            "L": -3,
            "N": -2,
            "Q": -1,
            "P": 7,
            "S": -1,
            "R": -2,
            "T": -1,
            "W": -4,
            "V": -2,
            "Y": -3,
            "X": -2,
            "Z": -1
        },
        "S": {
            "-": -4,
            "A": 1,
            "C": -1,
            "B": 0,
            "E": 0,
            "D": 0,
            "G": 0,
            "F": -2,
            "I": -2,
            "H": -1,
            "K": 0,
            "M": -1,
            "L": -2,
            "N": 1,
            "Q": 0,
            "P": -1,
            "S": 4,
            "R": -1,
            "T": 1,
            "W": -3,
            "V": -2,
            "Y": -2,
            "X": 0,
            "Z": 0
        },
        "R": {
            "-": -4,
            "A": -1,
            "C": -3,
            "B": -1,
            "E": 0,
            "D": -2,
            "G": -2,
            "F": -3,
            "I": -3,
            "H": 0,
            "K": 2,
            "M": -1,
            "L": -2,
            "N": 0,
            "Q": 1,
            "P": -2,
            "S": -1,
            "R": 5,
            "T": -1,
            "W": -3,
            "V": -3,
            "Y": -2,
            "X": -1,
            "Z": 0
        },
        "T": {
            "-": -4,
            "A": 0,
            "C": -1,
            "B": -1,
            "E": -1,
            "D": -1,
            "G": -2,
            "F": -2,
            "I": -1,
            "H": -2,
            "K": -1,
            "M": -1,
            "L": -1,
            "N": 0,
            "Q": -1,
            "P": -1,
            "S": 1,
            "R": -1,
            "T": 5,
            "W": -2,
            "V": 0,
            "Y": -2,
            "X": 0,
            "Z": -1
        },
        "W": {
            "-": -4,
            "A": -3,
            "C": -2,
            "B": -4,
            "E": -3,
            "D": -4,
            "G": -2,
            "F": 1,
            "I": -3,
            "H": -2,
            "K": -3,
            "M": -1,
            "L": -2,
            "N": -4,
            "Q": -2,
            "P": -4,
            "S": -3,
            "R": -3,
            "T": -2,
            "W": 11,
            "V": -3,
            "Y": 2,
            "X": -2,
            "Z": -3
        },
        "V": {
            "-": -4,
            "A": 0,
            "C": -1,
            "B": -3,
            "E": -2,
            "D": -3,
            "G": -3,
            "F": -1,
            "I": 3,
            "H": -3,
            "K": -2,
            "M": 1,
            "L": 1,
            "N": -3,
            "Q": -2,
            "P": -2,
            "S": -2,
            "R": -3,
            "T": 0,
            "W": -3,
            "V": 4,
            "Y": -1,
            "X": -1,
            "Z": -2
        },
        "Y": {
            "-": -4,
            "A": -2,
            "C": -2,
            "B": -3,
            "E": -2,
            "D": -3,
            "G": -3,
            "F": 3,
            "I": -1,
            "H": 2,
            "K": -2,
            "M": -1,
            "L": -1,
            "N": -2,
            "Q": -1,
            "P": -3,
            "S": -2,
            "R": -2,
            "T": -2,
            "W": 2,
            "V": -1,
            "Y": 7,
            "X": -1,
            "Z": -2
        },
        "X": {
            "-": -4,
            "A": 0,
            "C": -2,
            "B": -1,
            "E": -1,
            "D": -1,
            "G": -1,
            "F": -1,
            "I": -1,
            "H": -1,
            "K": -1,
            "M": -1,
            "L": -1,
            "N": -1,
            "Q": -1,
            "P": -2,
            "S": 0,
            "R": -1,
            "T": 0,
            "W": -2,
            "V": -1,
            "Y": -1,
            "X": -1,
            "Z": -1
        },
        "Z": {
            "-": -4,
            "A": -1,
            "C": -3,
            "B": 1,
            "E": 4,
            "D": 1,
            "G": -2,
            "F": -3,
            "I": -3,
            "H": 0,
            "K": 1,
            "M": -1,
            "L": -3,
            "N": 0,
            "Q": 3,
            "P": -1,
            "S": 0,
            "R": 0,
            "T": -1,
            "W": -3,
            "V": -2,
            "Y": -2,
            "X": -1,
            "Z": 4
        }
    }
};

function fact(n) {
    return n > 1 ? fact(n - 1) * n : 1
}

function div(val, by) {
    return (val - val % by) / by;
}

function subs(str) {
    let needle;
    [str, needle] = str.split('\n');
    str = str.trim();
    needle = needle.trim();
    let result = [];

    let ind = 0;

    while ((ind = str.indexOf(needle, ind)) > 0) {
        ind++;
        result.push(ind);
    }

    console.log(result.join(' '))
}

function hamm(str1, str2) {
    let count = 0;

    for (let i = 0; i < str1.length; i++) {
        count += str1[i] !== str2[i];
    }

    return count;
}

function swap(a, i, j) {
    let s = a[i];
    a[i] = a[j];
    a[j] = s;
}

function NextSet(a, n) {
    let j = n - 2;
    while (j != -1 && a[j] >= a[j + 1]) j--;
    if (j == -1)
        return false; // больше перестановок нет
    let k = n - 1;
    while (a[j] >= a[k]) k--;
    swap(a, j, k);
    let l = j + 1, r = n - 1; // сортируем оставшуюся часть последовательности
    while (l < r)
        swap(a, l++, r--);
    return true;
}

function Print(a) {
    console.log(a.join(' '))
}

function perm(n) {
    console.log(fact(n))
    let a = [];
    for (let i = 0; i < n; i++) {
        a[i] = i + 1;
    }
    Print(a, n);
    while (NextSet(a, n))
        Print(a, n);
}

function rev_comp(str) {
    return str.split('').reverse().join('').replace(/A|G|C|T/g, (m) => useful.DNAToDNA[m])
}

function orf(str) {
    let res = [];
    str = str.split('\n')[1].trim();

    res = res.concat(find_prot_dna(str));
    res = res.concat(find_prot_dna(str.substring(1)));
    res = res.concat(find_prot_dna(str.substring(2)));
    let compl = str.split('').reverse().join('').replace(/A|G|C|T/g, (m) => useful.DNAToDNA[m])


    res = res.concat(find_prot_dna(compl));
    res = res.concat(find_prot_dna(compl.substring(1)));
    res = res.concat(find_prot_dna(compl.substring(2)));
    console.log([...new Set(res)].join('\n'))

}

function matchOverlap(input, re) {
    var r = [], m;
    // prevent infinite loops
    if (!re.global) re = new RegExp(
        re.source, (re+'').split('/').pop() + 'g'
    );
    while (m = re.exec(input)) {
        re.lastIndex -= m[0].length - 1;
        r.push(m[1]);
    }
    return r;
}

function find_prot_dna(str) {
    let codons = str.match(/.{1,3}/g);
    let amins = codons.map(c => useful.codonToAminoDNA[c]).join('');
    let pat = /(M.*?)\*/g
    let results = matchOverlap(amins, pat)
    return results
}

function revp(str) {
    let res = [];
    str = str.split('\n')[1].trim();
    for (let i = 0; i < str.length; i++) {
        for (let j = 2; j <= 6; j++) {
            let sub = str.substr(i, j);
            let comp_sub = rev_comp(sub);
            if (str.substr(i + j, j) == comp_sub) {
                res.push(i + 1 + ' ' + j * 2)
            }
        }
    }
    console.log(res.join('\n'));
}

function grph(str) {
    let result = [];
    str = str.split('>').filter((e) => e);
    let data = str.map(row => {
        let d = row.substr(row.indexOf('\n')).replace('\n', '').trim();
        return {name: row.substring(0, row.indexOf('\n')), data: d}
    });
    data.forEach(elem1 => {
        data.forEach(elem2 => {
            if (elem1.name !== elem2.name) {
                if (elem1.data.substr(-3) == elem2.data.substr(0, 3)) {
                    result.push(elem1.name + ' ' + elem2.name);
                }
            }
        })
    });
    console.log(result.join('\n'))
}


function ba3b(str) {
    let [result, ...strings] = str.split('\n').filter(s => s.length).map(s => s.trim());

    strings.forEach(s => result +=s[s.length - 1]);
    console.log(result);
}

let arr;
let acc = '';
function long(str) {
    let data = str.split('>').filter(s => s).map(s => {
        [_, ...s] = s.split('\n').map(e => e.trim());
        return s.join('');
    })

    console.log(data);
    arr = data;
    console.log(find_overlaps());
    console.log(acc);
}

function  find_overlaps() {
    console.log('------------------------',arr, acc);
    if (arr.length == 0) {
        return acc;
    } else if (acc.length ==0) {
        acc = arr.shift();
        return find_overlaps(acc);
    } else {
        for (let i = 0; i < arr.length; i++) {
            let a = arr[i];
            let len = a.length;
            console.log(len, a);
            for (let p = len - 1; p >= 0; p--) {

                let q = len - p;
                let suff = a.substr(0, p)
                let pref = a.slice(q);
                console.log(p, q, acc, pref, acc.startsWith(pref), suff, acc.endsWith(suff), a, arr);
                // console.log(pref, suff);
                // console.log(acc.startsWith(pref), acc.endsWith(suff));
                if (acc.endsWith(suff) && suff) {
                    arr.splice(i, 1);
                    acc = acc + a.substr(p);
                    return find_overlaps()
                }
                if (acc.startsWith(pref) && pref) {
                    arr.splice(i, 1)
                    acc = a.substr(0, q) + acc;
                    return find_overlaps()
                }
                // ATTAGACCTGCCGGAAGACCTGCCGGAATAC
                // ATTAGACCTGCCGGAATAC
                // console.log(p, q, acc,
                //     a.substr(0, p), acc.endsWith(a.substr(0, p)),
                //     a.slice(q), acc.endsWith(a.slice(q)))
                // if (a.substr(0, p) && acc.startsWith(a.substr(0, p))) {
                //     // console.log('-----------',a.substr(p), acc)
                //     arr.splice(i, 1);
                //     acc = a.substr(p) + acc;
                //     return find_overlaps()
                // }
                // if (a.slice(q) && acc.endsWith(a.slice(q))) {
                //     // console.log('-----------',a.slice(q), acc)
                //     arr.splice(i, 1);
                //     acc = acc + a.slice(q)
                //     return find_overlaps();
                // }
            }
        }
    }
}

function long(str) {
    let data = str.split('>').filter(s=>s).map(s=>{
        [_, ...s] = s.split('\n').map(e=> e.trim());
        return s.join('');
    });
    data.sort((s1, s2)=> s2.length - s1.length);
    // console.log(data);


    let result = data.shift();
    let count = data.length;
    while (data.length > 0 && count > 0) {
        let len_ar = data.map(s => s.length);
        // console.log(len_ar);
        loop1:
            for (var i = len_ar.shift() - 1; i > 0; i--) {
                // console.log(i);
                let pref_ar = data.map(e => e.substring(e.length-i))
                let suff_ar = data.map(e => e.substring(0, i))
                // data.forEach(e => console.log(e.length-i, e.substring(e.length-i), result, e.substring(0, i)))
                for (let j = 0; j < pref_ar.length; j++) {
                    // let pref = pref_ar[j];
                    let pref = data[j]
                    let suff = suff_ar[j];
                    if (result.startsWith(pref)) {
                        result = data[j].substring(0, data[j].length - i) + result;
                        data.splice(i, 1);
                        count--;
                        break loop1;
                    }
                    if (result.endsWith(suff)) {
                        result = result + data[j].substring(i);
                        data.splice(i, 1);
                        count --;
                        break loop1;
                    }
                }

            }
    }
    console.log(result)
}


function corr(str) {
    let data = str.split('>').filter(s=>s).map(s=>{
        [_, ...s] = s.split('\n').map(e=> e.trim());
        return s.join('');
    });
    let data_copy = data.slice();
    let correct = [];

    while (data.length > 0) {
        let s = data.shift();
        let s_rev = rev_comp(s);
        for (let i = 0; i < data.length; i++) {
            let d = hamm(s, data[i]);
            let d_rev = hamm(s_rev, data[i]);
            if (Math.min(d, d_rev) === 0) {
                correct.push(s)
            }
        }
    }
    correct = [...new Set(correct)];
    data_copy.forEach(s => {
        let s_rev = rev_comp(s);
        for (let i = 0; i < correct.length; i++) {
            let d = hamm(s, correct[i]);
            let d_rev = hamm(s_rev, correct[i]);
            if (d === 1) {
                console.log(s + '->' + correct[i]);
                return;
            }
            if (d_rev === 1) {
                console.log(s + '->' + rev_comp(correct[i]))
            }
        }
    })
}

function ba3d(str) {
    let [n, text] = str.split('\n').map(s => s.trim());
    n = parseInt(n);
    let v = [];
    for (let i = 0; i < str.length - n; i++) {
        let s = text.substr(i, n);
        s.length === n  ? v.push([s.substr(0, n-1), s.substr(1)]) : '';
    }
    // console.log(v)
    while (v.length > 0) {
        let [n1, n2] = v.shift();
        let left = [n2];
        let i = 0;
        while (i < v.length) {
            if (n1 === v[i][0]) {
                left.push(v[i][1]);
                v.splice(i, 1)
            } else {
                i++;
            }
        }
        console.log(n1, '->', left.join(','))
    }
}

function ba3e(str) {
    let v = str.split('\n').filter(e=>e).map(e => {
        e = e.trim();
        return [e.substr(0, e.length - 1), e.substr(1)]
    });
    while (v.length > 0) {
        let [n1, n2] = v.shift();
        let left = [n2];
        let i = 0;
        while (i < v.length) {
            if (n1 === v[i][0]) {
                left.push(v[i][1]);
                v.splice(i, 1)
            } else {
                i++;
            }
        }
        console.log(n1, '->', left.join(','))
    }

}

function ba3f(str) {
    let rebra = str.split('\n').filter(e=>e).reduce((res, row) => {
        let [a, bs] = row.split('->');
        a = a.trim();
        bs = bs.split(',').map(b=>b.trim());
        bs.map(b=> res.push([a, b]))
        // res[a] = bs;
        return res;
    }, []);
    // console.log(rebra);
    let result = [];
    console.log(find_euler(rebra).reverse().join(' -> '))

}

function find_euler(graph) {
    let stack = [];
    let tour = [];

    stack.push(graph[0][0]);

    while (stack.length > 0) {
        // console.log(stack, graph, tour)
        // console.log()
        let v = stack[stack.length - 1];


        let d = get_degr(v, graph);

        if (d === 0) {
            stack.pop();
            tour.push(v)
        } else {
            let [index, edge] = get_ind_edge(v, graph);
            // console.log(v, d, index, edge)
            graph.splice(index, 1);
            stack.push(v === edge[0] ? edge[1] : edge[0])
        }

    }
    return tour;
}

function get_degr(v, graph) {
    let d = 0;
    for (let [x, y] of graph) {
        if (x === v) {
            d += 1
        }
    }
    return d;
}

function get_ind_edge(v, graph) {

    for (let i = 0; i < graph.length; i++) {
        if (v === graph[i][0] || v === graph[i][1]) {
            return [i, graph[i]]
        }
    }
}

function levenshtein(s1, s2, costs) {
    let i, j, l1, l2, flip, ch, chl, ii, ii2, cost, cutHalf;
    l1 = s1.length;
    l2 = s2.length;

    costs = costs || {};
    let cr = costs.replace || 1;
    let cri = costs.replaceCase || costs.replace || 1;
    let ci = costs.insert || 1;
    let cd = costs.remove || 1;

    cutHalf = flip = Math.max(l1, l2);

    let minCost = Math.min(cd, ci, cr);
    let minD = Math.max(minCost, (l1 - l2) * cd);
    let minI = Math.max(minCost, (l2 - l1) * ci);
    let buf = new Array((cutHalf * 2) - 1);

    for (i = 0; i <= l2; ++i) {
        buf[i] = i * minD;
    }

    for (i = 0; i < l1; ++i, flip = cutHalf - flip) {
        ch = s1[i];
        chl = ch.toLowerCase();

        buf[flip] = (i + 1) * minI;

        ii = flip;
        ii2 = cutHalf - flip;

        for (j = 0; j < l2; ++j, ++ii, ++ii2) {
            cost = (ch === s2[j] ? 0 : (chl === s2[j].toLowerCase()) ? cri : cr);
            buf[ii + 1] = Math.min(buf[ii2 + 1] + cd, buf[ii] + ci, buf[ii2] + cost);
        }
    }
    return buf[l2 + cutHalf - flip];
}

function ba5g(str) {
    let [s1, s2] = str.split('\n').map(e => e.trim());
    console.log(levenshtein(s1, s2))
}

function ba5k(str) {
    let [s1, s2] = str.split('\n').map(e => e.trim());
    console.log(middleColumn(s1, s2))
}

function middleColumn(v, w, indel) {
    let S = [[],[]];
    for(let i = 0; i <= v.length; ++i) {
        S[0][i] = i * -1 * indel;
    }
    S[1][0] = indel;
    let backtrack = [];
    for(var j = 1; j <= w.length/2; ++j) {
        for(var i = 0; i <= v.length; ++i) {
            if(i == 0) {
                S[1][i] = -1*j*indel;
            }
            else {
                let opt1 = S[0][i-1] + useful.blosum62[v[i-1]][w[j-1]];
                let opt2 = S[0][i] - indel;
                let opt3 = S[1][i-1] - indel;
                S[1][i] = opt1;
                backtrack[i] = 0;
                if(opt2 > S[1][i]) {
                    S[1][i] = opt2;
                    backtrack[i] = 1;
                }
                if(opt3 > S[1][i]) {
                    S[1][i] = opt3;
                    backtrack[i] = 2;
                }
            }
        }
        if(j != w.length/2) {
            S[0] = S[1];
            S[1] = [];
        }
    }
    console.log(S``,backtrack);
}

function ba5e(str) {
    let [s1, s2] = str.split('\n').map(e => e.trim());
    allign_prot_go_ge(s1, s2);
}

function allign_prot(s,  t, indel = 5) {
    let S = [[]];
    let opt = [[]]; // 1 = right, 2 = down, 3 = diag
    S[0][0] = 0;
    for(let i = 1; i <= s.length; ++i) {
        S[i] = [];
        opt[i] = [];
        S[i][0] = S[i-1][0] - indel;
        opt[i][0] = 1;
    }
    for(let j = 1; j <= t.length; ++j) {
        S[0][j] = S[0][j-1] - indel;
        opt[0][j] = 2;
    }
    for(let i = 1; i <= s.length; ++i) {
        for(let j = 1; j <= t.length; ++j) {
            let opt1 = S[i][j-1] - indel; // right
            let opt2 = S[i-1][j] - indel; // down
            let opt3 = S[i-1][j-1] + useful.blosum62[s[i-1]][t[j-1]]; // diag
            S[i][j] = opt1;
            opt[i][j] = 1;
            if(opt2 > S[i][j]) {
                S[i][j] = opt2;
                opt[i][j] = 2;
            }
            if(opt3 > S[i][j]) {
                S[i][j] = opt3;
                opt[i][j] = 3;
            }
        }
    }
    let out = ['', ''];
    let count = S[s.length][t.length];

    let i = s.length;
    let j = t.length;
    while(i > 0 && j > 0) {
        if(opt[i][j] == 1) { // right
            out[0] = '-' + out[0];
            out[1] = t[j-- - 1] + out[1];
        }
        else if(opt[i][j] == 2) { // down
            out[0] = s[i-- - 1] + out[0];
            out[1] = '-' + out[1];
        }
        else if(opt[i][j] == 3) { // diag
            out[0] = s[i-- - 1] + out[0];
            out[1] = t[j-- - 1] + out[1];
        }
    }
    if(i > 0) {
        out[0] = s.substring(0,i) + out[0];
        let add = "";
        for(let x = 0; x < i; ++x) {
            add += '-';
        }
        out[1] = add + out[1];
    }
    if(j > 0) {
        out[1] = t.substring(0,j) + out[1];
        let add = "";
        for(let x = 0; x < j; ++x) {
            add += '-';
        }
        out[0] = add + out[0];
    }
    console.log(count)
    console.log(out[0])
    console.log(out[1])
    return out[0];
}

function allign_prot_go_ge(s, t, go = 11, ge = 1) {
    // Uses BLOSUM62. PASS go AND ge AS POSITIVE INTEGER, NOT NEGATIVE!
    let S = [[[]],[[]],[[]]]; // S[0] = lower, S[1] = middle, S[2] = upper
    let opt = [[[]],[[]],[[]]]; // 1 = right, 2 = down, 3 = diag
    for(let i = 1; i <= s.length; ++i) {
        S[0][i] = [];
        S[1][i] = [];
        S[2][i] = [];
        opt[0][i] = [];
        opt[1][i] = [];
        opt[2][i] = [];
        S[0][i][0] =  (go + (i-1)*ge);
        S[1][i][0] =  (go + (i-1)*ge);
        S[2][i][0] = go;
        opt[0][i][0] = 0;
        opt[1][i][0] = 0;
        opt[2][i][0] = 0;
    }
    for(let j = 1; j <= t.length; ++j) {
        S[0][0][j] = go;
        S[1][0][j] =  (go + (j-1)*ge);
        S[2][0][j] =  (go + (j-1)*ge);
        opt[0][0][j] = 1;
        opt[1][0][j] = 1;
        opt[2][0][j] = 1;
    }
    for(let i = 1; i <= s.length; ++i) {
        for(let j = 1; j <= t.length; ++j) {
            let low1 = S[0][i-1][j] - ge;
            let low2 = S[1][i-1][j] - go;
            if(low1 > low2) {
                S[0][i][j] = low1;
                opt[0][i][j] = 0;
            }
            else {
                S[0][i][j] = low2;
                opt[0][i][j] = 1;
            }

            let up1 = S[2][i][j-1] - ge;
            let up2 = S[1][i][j-1] - go;
            if(up1 > up2) {
                S[2][i][j] = up1;
                opt[2][i][j] = 0;
            }
            else {
                S[2][i][j] = up2;
                opt[2][i][j] = 1;
            }

            let opt1 = S[0][i][j];
            let opt2 = S[1][i-1][j-1] + useful.blosum62[s[i-1]][t[j-1]];
            let opt3 = S[2][i][j];
            S[1][i][j] = opt1;
            opt[1][i][j] = 0;
            if(opt2 > S[1][i][j]) {
                S[1][i][j] = opt2;
                opt[1][i][j] = 1;
            }
            if(opt3 > S[1][i][j]) {
                S[1][i][j] = opt3;
                opt[1][i][j] = 2;
            }
        }
    }
    let i = s.length;
    let j = t.length;
    let bestSIJ = 0;
    let best = S[0][i][j];
    if(S[1][i][j] > best) {
        best = S[1][i][j];
        bestSIJ = 1;
    }
    if(S[2][i][j] > best) {
        best = S[2][i][j];
        bestSIJ = 2;
    }
    let out = [];
    out[0] = "" + best;
    out[1] = "";
    out[2] = "";
    while(i > 0 && j > 0) {
        if(bestSIJ == 0) {
            if(opt[0][i][j] == 1) {
                bestSIJ = 1;
            }
            out[1] = s[i-- - 1] + out[1];
            out[2] = '-' + out[2];
        }
        else if(bestSIJ == 1) {
            if(opt[1][i][j] == 0) {
                bestSIJ = 0;
            }
            else if(opt[1][i][j] == 2) {
                bestSIJ = 2;
            }
            else {
                out[1] = s[i-- - 1] + out[1];
                out[2] = t[j-- - 1] + out[2];
            }
        }
        else {
            if(opt[2][i][j] == 1) {
                bestSIJ = 1;
            }
            out[1] = '-' + out[1];
            out[2] = t[j-- - 1] + out[2];
        }
    }
    if(i > 0) {
        out[1] = s.substring(0,i) + out[1];
        let add = "";
        for(let x = 0; x < i; ++x) {
            add += '-';
        }
        out[2] = add + out[2];
    }
    if(j > 0) {
        out[2] = t.substring(0,j) + out[2];
        let add = "";
        for(let x = 0; x < j; ++x) {
            add += '-';
        }
        out[1] = add + out[1];
    }
    console.log(out[0]);
    console.log(out[1]);
    console.log(out[2]);
    return out;
}


fs.readFile('tmp.txt', (err, data) => {
// fs.readFile('rosalind_ba5j.txt', (err, data) => {
    ba5k(data.toString());
});
