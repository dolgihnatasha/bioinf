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
}

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

    // let [str1, str2] = str.split('\n');
    // str1 = str1.trim();
    // str2 = str2.trim();
    let count = 0;
    // let count_rev = 0;
    // let str2_rev = rev_comp(str2);

    for (let i = 0; i < str1.length; i++) {
        count += str1[i] !== str2[i];
        // count_rev += str1 !== str2_rev[i];
    }

    return count;
    // return min(count, count_rev);
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
    // console.log(data);
    data.forEach(elem1 => {
        data.forEach(elem2 => {
            if (elem1.name !== elem2.name) {
                if (elem1.data.substr(-3) == elem2.data.substr(0, 3)) {
                    result.push(elem1.name + ' ' + elem2.name);
                }
            }
        })
    })
    console.log(result.join('\n'))
}


function ba3b(str) {
    let [result, ...strings] = str.split('\n').filter(s => s.length).map(s => s.trim());

    strings.forEach(s => result +=s[s.length - 1]);
    console.log(result);
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

// while (rebra.length > 0) {
//     let v = result[result.length - 1];
//
//     loop:
//         for (let i = 0; i < rebra.length; i++) {
//             console.log(v, rebra[i][0], result)
//             if (rebra[i][0] == v) {
//                 result.push(rebra[i][1]);
//                 rebra.splice(i, 1);
//                 break loop;
//             }
//         }
// }
// function findEulerPath(v) {
//     if (rebra[v]) {
//
//         // console.log(v)
//         let u = rebra[v].slice();
//         delete rebra[v];
//         u.forEach(findEulerPath);
//     }
//     // console.log(v, result.length, result, JSON.stringify(rebra))
//     result.push(v)
//
//     return;
// }
// findEulerPath(0)
//
// console.log(result.reverse().join(' -> '))


// fs.readFile('tmp.txt', (err, data) => {
fs.readFile('rosalind_ba3f.txt', (err, data) => {
    ba3f(data.toString());
});
