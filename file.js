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

function hamm(str) {

    let [str1, str2] = str.split('\n');
    str1 = str1.trim();
    str2 = str2.trim();
    let count = 0;

    for (let i = 0; i < str1.length; i++) {
        count += str1[i] !== str2[i];
    }

    console.log(count)
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


fs.readFile('tmp.txt', (err, data) => {
    long(data.toString());
});
