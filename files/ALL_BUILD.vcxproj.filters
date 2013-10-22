from django.db import models
from multiprocessing import Process

from sufpref import sf
import string

import json
import base64

fields = [
    "title", 
    "tag", 
    "text_length", 
    "no", 
    "alphabet", 
    "order",
    "order_value", 
    "order_minus2_1", 
    "order_minus1_1", 
    "order_minus1_2", 
    "order_0_1", 
    "order_gt0_1",
    "mark_ini_probs",
    "pattern",
    "wordslist", 
    "nwords", 
    "wordlen", 
    "bernprob", 
    "cutoff", 
    "pssm1", 
    "pssm2", 
    "footprints", 
    "word", 
    "nreplace", 
    "constpositions", 
    "cons_alp", 
    "cons_word"
]

class Job(models.Model):
    ##### common fields
    # primary key
    id = models.AutoField(primary_key=True)
    # starter IP
    ip = models.IPAddressField()
    # start time
    started = models.DateTimeField(auto_now=False, auto_now_add=True)
    # end time
    finished = models.DateTimeField(auto_now=True, auto_now_add=False)
    # done
    done = models.BooleanField(default=False)
    ##### comment fields
    # job title
    title = models.TextField()
    # username tag
    tag = models.TextField()
    ##### specific fields
    text_length = models.IntegerField()
    no = models.IntegerField()
    alphabet = models.TextField()
    # order-dependend values
    order = models.IntegerField()
    hhm_probs = models.TextField() # Json-encoded 3D array
    dhhm_probs = models.TextField() # Json-encoded 2D array
    dhhm_trans = models.TextField() # Json-encoded 2D array
    bern_probs = models.TextField() # Json-encoded 1D array
    mark_probs = models.TextField() # Json-encoded 2D array
    mark_ini_probs = models.TextField() # Json-encoded 1D array
    # pattern-depended values
    pattern = models.IntegerField()
    wordslist = models.TextField() # Json-encoded string list
    nwords = models.IntegerField()
    wordlen = models.IntegerField()
    bernprob = models.TextField() # Json-encoded 1D array
    cutoff = models.FloatField()
    pssm1 = models.TextField() # Json-encoded 2D array
    pssm2 = models.TextField() # Json-encoded 2D array
    footprints = models.TextField() # Json-encoded string list
    word = models.TextField()
    nreplace = models.IntegerField()
    constpositions = models.TextField() # Json-encoded 1D array
    cons_alp = models.TextField() # Json-encoded dict: char -> string
    cons_word = models.TextField()
    result = models.FloatField()
    error = models.IntegerField()
    report = models.TextField()
    res_words = models.TextField() # Json-encoded string list


def parse_2D_matrix(text, is_integer):
    matrix = []
    error = False
    lines = text.split("\n")
    for line in lines:
        if error:
            break
        if len(line.strip())==0:
            continue
        terms = line.split()
        row = []
        for term in terms:
            if error:
                break
            try:
                if is_integer:
                    value = int(term)
                else:
                    value = float(term)
            except:
                value = 0.0
                error = True
            row.append(value)
        matrix.append(row)
    return matrix, error

def create_db_entry(args):
    job = Job()
    job.title = args["title"].strip()
    job.tag = args["tag"].strip()
    job.ip = args["ip"]
    errors = []
    try:
        job.text_length = int(args["text_length"])
    except:
        job.text_length = -1
        if len(args["text_length"])==0:
            errors += ["Text length is empty"]
        else:
            errors += ["Text length is not integer value"]
    try:
        job.no = int(args["no"])
    except:
        job.no = -1
        if len(args["no"])==0:
            errors += ["Number of occurences is empty"]
        else:
            errors += ["Number of occurences is not integer value"]
    job.alphabet = args["alphabet"].strip().replace(" ","").replace(",","").replace(" ","").upper()
    if len(job.alphabet)==0:
        errors += ["Alphabet is empty"]
    try:
        job.order = int(args["order"])
        if job.order==1:
            job.order = int(args["order_value"])
        
    except:
        job.order = -10;
        errors += ["Order value in not integer"]
    hhm_lines = args["order_minus2_1"].split("\n")
    hhm = []
    hhm_matrix = []
    hhm_error = False
    for line in hhm_lines:
        if hhm_error:
            break
        if len(line.strip())==0:
            continue
        terms = line.split()
        if len(terms)==1 and terms[0][0]=='-':
            if len(hhm_matrix)>0:
                hhm.append(hhm_matrix)
                hhm_matrix = []
            continue
        row = []
        for term in terms:
            if hhm_error:
                break
            try:
                value = float(term)
            except:
                value = 0.0
                hhm_error = True
                if job.order==-2:
                    errors += ["Probability distribution of HHM contains non floating-point value"]
            row.append(value)
        hhm_matrix.append(row)
    if len(hhm_matrix)>0:
        hhm.append(hhm_matrix)
    
    if (job.order==-2):
        # Check:
        # 1) matrices are the same size
        # 2) rows count = matrices count
        # 3) columns count = alphabet size
        matrix_no = 0
        matrices_count = len(hhm)
        for hhm_matrix in hhm:
            matrix_no += 1
            rows_count = len(hhm_matrix)
            if rows_count!=matrices_count:
                errors += ["HHM Probs: matrix "+str(matrix_no)+" rows count doesn't match matrices count"]
            else:
                for row in hhm_matrix:
                    if len(row)!=len(job.alphabet) and len(job.alphabet)>0:
                        errors += ["HHM Probs: matrix "+str(matrix_no)+" columns count doesn't match alphabet size"]
                        break
    
    job.hhm_probs = json.dumps(hhm)
    
    dhhm_p, error = parse_2D_matrix(args["order_minus1_1"], False)
    job.dhhm_probs = json.dumps(dhhm_p)
    if job.order==-1 and error:
        errors += ["Probabilities of Determenistic HHM contains non floating-point value"]
        
    dhhm_t, error = parse_2D_matrix(args["order_minus1_2"], True)
    job.dhhm_trans = json.dumps(dhhm_t)
    if job.order==-1 and error:
        errors += ["Transitions of Determenistic HHM contains non integer value"]
        
    if job.order==-1:
        # Check:
        # 1) dhhm_p and dhhm_t sizes are equal
        # 2) dhhm_p/dhhm_t columns count == alphabet size
        if len(dhhm_p)!=len(dhhm_t):
            errors += ["Determenistic HHM: Probabilities and transitions matrices have different rows count"]
        else:
            for row in dhhm_p:
                if len(row)!=len(job.alphabet) and len(job.alphabet)>0:
                    errors += ["Determenistic HHM: Probabilities matrix columns count doesn't match alphabet size"]
                    break
            for row in dhhm_t:
                if len(row)!=len(job.alphabet) and len(job.alphabet)>0:
                    errors += ["Determenistic HHM: Transitions matrix columns count doesn't match alphabet size"]
                    break
    
    bern_prob_line = args["order_0_1"].strip()
    terms = bern_prob_line.split()
    bern_probs = []
    bern_probs_error = False
    for term in terms:
        if bern_probs_error:
            break
        try:
            value = float(term)
        except:
            value = 0.0
            bern_probs_error = True
            if job.order==0:
                errors += ["Probability distribution of Bernoulli contains non floating-point value"]
        bern_probs.append(value)
    
    job.bern_probs = json.dumps(bern_probs)
    
    if job.order==0:
        # Check
        # 1) columns count == alphabet size
        if len(bern_probs)!=len(job.alphabet) and len(job.alphabet)>0:
            errors += ["Bernoulli: Probability distribution elements count doesn't match alphabet size"]
    
    mark_probs, error = parse_2D_matrix(args["order_gt0_1"], False)    
    job.mark_probs = json.dumps(mark_probs)
    
    if (job.order>0 or job.order==-10) and error:
        errors += ["Probability distribution of Markovian model contains non floating-point value"]
        
    mark_ini_probs_line = args["mark_ini_probs"].strip()
    terms = mark_ini_probs_line.split()
    mark_ini_probs = []
    mark_ini_probs_error = False
    for term in terms:
        if mark_ini_probs_error:
            break
        try:
            value = float(term)
        except:
            value = 0.0
            mark_ini_probs_error = True
            if job.order>0 or job.order==-10:
                errors += ["Initial probabilities of Markovian model contains non floating-point value"]
        mark_ini_probs.append(value)
        
    job.mark_ini_probs = json.dumps(mark_ini_probs)
            
    if job.order>0 or job.order==-10:
        # Check
        # 1) columns count == alphabet size
        # 2) rows count == alphabet size ** order value
        # 3) mark_ini_probs (if not empty) size = AlpSize**order
        if job.order>0 and len(job.alphabet)>0:
            must_be_rows_count = len(job.alphabet)**job.order
            if len(mark_probs)!=must_be_rows_count:
                errors += ["Markovian Model: probability distribution matrix rows count doesn't equal alphabet size in power of order value"]
            else:
                for row in mark_probs:
                    if len(row)!=len(job.alphabet) and len(job.alphabet)>0:
                        errors += ["Markovian Model: probability distribution matrix columns count doesn't match alphabet size"]
                        break
        if job.order>0 and len(job.alphabet)>0 and len(mark_ini_probs)>0:
            must_be = len(job.alphabet)**job.order
            if len(mark_ini_probs)!=must_be:
                errors += ["Markovian Model: initial probabilities count doesn't equal alphabet size in power of order value"]
    
    try:
        job.pattern = int(args["pattern"])
    except:
        job.pattern = -10
    
    wordslist = filter(lambda x: len(x.strip())>0, args["wordslist"].split("\n"))
    wordslist = map(lambda x: x.replace(" ","").replace("\n","").replace("\r","").upper(), wordslist)
    job.wordslist = json.dumps(wordslist)
    
    if job.pattern==0:
        # Check:
        # 1) forall word in list: |word|==constant
        # 2) forall word in list: forall char in word: char in alphabet
        word_size = -1
        for word in wordslist:
            if word_size==-1:
                word_size = len(word)
            else:
                if len(word)!=word_size:
                    errors += ["List of words: different word sizes"]
                    break
        word_no = 0
        for word in wordslist:
            word_no += 1
            for ch in word:
                if not ch in job.alphabet:
                    errors += ["List of words: word in line "+str(word_no)+" contains symbol out of alphabet"]
                    break
    
    try:
        job.nwords = int(args["nwords"])
    except:
        job.nwords = 0
        if job.pattern==1:
            errors += ["Random generated pattern: number of words is not integer"]
    
    try:
        job.wordlen = int(args["wordlen"])
    except:
        job.wordlen = 0
        if job.pattern==1:
            errors += ["Random generated pattern: length of word is not integer"]
    
    bernprobs = args["bernprob"].split()
    bp = []
    for p in bernprobs:
        try:
            value = float(p)
        except:
            errors += ["Random generated pattern: not floating-point value in list"]
            value = 0.0
        bp.append(value)
    job.bernprob = json.dumps(bp)
    
    if job.pattern==1:
        if len(bernprobs)!=len(job.alphabet) and len(job.alphabet)>0:
            errors += ["Random generated pattern: Probability distribution elements count doesn't match alphabet size"]
    
    try:
        job.cutoff = float(args["cutoff"])
    except:
        job.cutoff = 0.0
        if job.pattern==2:
            errors += ["PSSM and cut-off: cut-off is not floating-point value"]
    
    pssm1,error = parse_2D_matrix(args["pssm1"], False)
    job.pssm1 = json.dumps(pssm1)
    
    if job.pattern==2 and error:
        errors += ["PSSM and cut-off: PSSM contains non floating-point value"]
        
    if job.pattern==2:
        # Check
        # 1) pssm columns count == alphabet size
        # 2) pssm rows count >= 1
        if len(pssm1)==0:
            errors += ["PSSM and cut-off: PSSM matrix is empty"]
        for row in pssm1:
            if len(row)!=len(job.alphabet) and len(job.alphabet)>0:
                errors += ["PSSM and cut-off: PSSM columns count doesn't match alphabet size"]
    
    pssm2,error = parse_2D_matrix(args["pssm2"], False)
    job.pssm2 = json.dumps(pssm2)
    
    if job.pattern==3 and error:
        errors += ["PSSM and footprints: PSSM contains non floating-point value"]
        
    if job.pattern==3:
        # Check
        # 1) pssm columns count == alphabet size
        if len(pssm2)==0:
            errors += ["PSSM and footprints: PSSM matrix is empty"]
        for row in pssm2:
            if len(row)!=len(job.alphabet) and len(job.alphabet)>0:
                errors += ["PSSM and footprints: PSSM columns count doesn't match alphabet size"]
    
    footprints = filter(lambda x: len(x.strip())>0, args["footprints"].split("\n"))
    footprints = map(lambda x: x.replace(" ","").upper(), footprints)
    job.footprints = json.dumps(footprints)
    
    if job.pattern==3 and len(job.alphabet)>0:
        word_no = 0
        for word in footprints:
            word_no +=1
            for ch in word:
                if not ch in job.alphabet:
                    errors += ["PSSM and footprints: footprint at line "+str(word_no)+" contains a symbol out of alphabet"]
                    break
    
    job.word = args["word"].strip().replace(" ","").upper()
    
    if job.pattern==4:
        for ch in job.word:
            if not ch in job.alphabet:
                errors += ["Word and mismatches: word contains a symbol out of alphabet"]
                break;
    
    try:
        job.nreplace = int(args["nreplace"])
    except:
        job.nreplace = 0
        if job.pattern==4:
            errors += ["Word and mismatches: maximum number of mismatches is not integer value"]
    
    constpositions_s = args["constpositions"].split()
    cps = []
    for cp in constpositions_s:
        cp = cp.strip()
        if len(cp)==0:
            continue
        try:
            value = int(cp)
        except:
            value = 0
            if job.pattern==4:
                errors += ["Word and mismatches: unchangeable positions contains non-integer value"]
            break
        cps.append(value)
    job.constpositions = json.dumps(cps)
    
    consalp_lines = args["cons_alp"].replace("[","{").replace("]","}").replace("(","{").replace(")","}").split("\n")
    consalp = dict()
    for line in consalp_lines:
        lexems = line.split("=")
        if (len(lexems))!=2:
            continue
        left = lexems[0].strip()
        right = lexems[1].replace(" ","").replace("{","").replace("}","")
        terms = right.split(",")
        terms = map(lambda x: x.strip(), terms)
        consalp[left] = terms
    
    if job.pattern==5:
        # 1. Check for correct symbols
        for key in consalp:
            if len(key)!=1:
                errors += ["Illegal symbol in left part of consensus alphabet: '"+key+"'"]
            right_part = consalp[key]
            for symbol in right_part:
                if len(symbol)!=1:
                    errors += ["Illegal symbol in right part of consensus alphabet: '"+symbol+"'"]
                else:
                    if len(job.alphabet)>0:
                        if not symbol in job.alphabet:
                            errors += ["Symbol in right part of consensus out of alphabet: '"+symbol+"'"]
        # 2. Check for varity of symbols
        for key in consalp:
            allsyms = string.join(consalp[key], "")
            badsyms = ""
            for sym in allsyms:
                if allsyms.count(sym)>1:
                    badsyms += sym
                    if not sym in badsyms:
                        errors += ["Duplicate symbol in right part of consensus out of alphabet: '"+sym+"'"]
            
    
    job.cons_alp = json.dumps(consalp)
    job.cons_word = args["cons_word"].strip()
    
    job.result = 0.0
    job.error = 0
    
    if len(errors)==0:
        job.save()
        id = job.id
    else:
        id = -1
    #assert False
    return id, errors
    
def run_sf_worker(job_id):
    job = Job.objects.get(id=job_id)
    local_error = sf.set_input_data(int(job.order),
                                  int(job.pattern),
                                  int(job.text_length),
                                  int(job.no))
    job.error = local_error
    
    
    if job.error!=0:
        job.done = True
        job.save()
        return
    
    alp_mas = str(job.alphabet)
    
    order = job.order
    
    if   order==-2:
        nd_hhm_probs = json.loads(job.hhm_probs)
        num_all_states = len(nd_hhm_probs)
        local_error = sf.analis_alp_hhm_data(alp_mas, num_all_states, nd_hhm_probs)
    elif order==-1:
        d_hhmprobs = json.loads(job.dhhm_probs)
        d_hhmtrans = json.loads(job.dhhm_trans)
        num_all_states = len(d_hhmprobs)
        local_error = sf.analis_alp_dhhm_data(alp_mas, num_all_states, d_hhmprobs, d_hhmtrans)
    elif order== 0:
        bern_prob = json.loads(job.bern_probs)
        local_error = sf.analis_alp_bern_data(alp_mas, bern_prob)
    elif order > 0:
        markov_probs = json.loads(job.mark_probs)
        mark_ini_probs = json.loads(job.mark_ini_probs)
        local_error = sf.analis_alp_mark_data(alp_mas, markov_probs, mark_ini_probs)
    job.error  = local_error
     
    
   
    if job.error!=0:
        job.done = True
        job.save()
        return
    
    pattern = job.pattern
    
    if   pattern==0:
        words_list = json.loads(job.wordslist)
        local_error = sf.analis_pattern_data_0(words_list)
    elif pattern==1:
        nwords = int(job.nwords)
        wordlen = int(job.wordlen)
        bern_prob = json.loads(job.bernprob)
        local_error = sf.analis_pattern_data_1(nwords, wordlen, bern_prob)
    elif pattern==2:
        pssm = json.loads(job.pssm1)
        thr = float(job.cutoff)
        wordlen = len(pssm)
        local_error = sf.analis_pattern_data_2_3(wordlen, [], pssm, thr)
    elif pattern==3:
        pssm = json.loads(job.pssm1)
        footprints = json.loads(job.footprints)
        wordlen = len(pssm)
        local_error = sf.analis_pattern_data_2_3(wordlen, footprints, pssm, 0.0)
    elif pattern==4:
        motif = str(job.word)
        nreplace = int(job.nreplace)
        constpositions = json.loads(job.constpositions)
        local_error = sf.analis_pattern_data_4(motif, nreplace, constpositions)
    elif pattern==5:
        consensus = str(job.cons_word)
        cons_alp_d = json.loads(job.cons_alp)
        cons_alp = []
        for key in cons_alp_d.keys():
            line = key+""+string.join(cons_alp_d[key], "")+""
            cons_alp.append(line)
        
        local_error = sf.analis_pattern_data_5(consensus, cons_alp)
    job.error = local_error
    
    
    
    if job.error!=0:
        job.done = True
        job.save()
        return
    
    
    local_error, job.result, local_report, local_res_words = sf.main()
    job.error = local_error
    job.report = local_report
    job.res_words = local_res_words
    
    
    
    job.done = True
    job.save()
    return
    
def run_sf(job_id):
    proc = Process(target=run_sf_worker, args=(job_id,))
    proc.start()
    #run_sf_worker(job_id)
    
ERRORS = [
    ""
    "Wrong number of parameters associated with the alphabet",
    "Error of opening of file with alphabet description	",
    "Invalid value of model order",
    "Incorrect distribution of letters. Frequences of letters must be real numbers",
    "Invalid distribution. Number of frequences must be equal to the size of alphabet",
    "Invalid distribution. Incorrect size of matrix with frequences of letters",
    "Incorrect alphabet distribution",
    "Invalid distribution. Incorrect size of matrix with transition states",
    "Wrong number of parameters associated with the pattern",
    "Incorrect parameters associated with the pattern",
    "Error of opening of file with pattern description",
    "Incorrect words in the pattern. Letters of words must be in the given alphabet",
    "Incorrect length of the words in the pattern",
    "Error of opening of file with PSSM",
    "Error of opening of file with footprints",
    "Error of opening of file with consensus alphabet description",
    "Incorrect description of consensus alphabet",
    "Incorrect letters of the consensus",
    "Wrong number of parameters associated with text length or number of occurences",
    "Incorrect parameters associated with text length or number of occurences",
    "Wrong number of parameters associated with the name of output file",
    "Wrong number of probabilities associated with Bernoulli probabilities",
    "Incorrect number of words in the pattern",
    "Incorrect length of the pattern",
    "Incorrect size of the alphabet",
    "Incorrect order of the probability model",
    "Incorrect number of states of the HHM",
    "Incorrect size of the alphabet",
    "The number of words must be smaller then  alphabet size in power of the pattern length",
    "Incorrect PSSM description",
    "Some error was caused on the previous stages",
    "Incorrect using of the function. This function can not be used with this probability model",
    "Incorrect using of the function. This function can non be used for such pattern description",
]
    
def get_status(job_id):
    job = Job.objects.get(id=job_id)
    result = ""
    if job.done:
        error = job.error
        if error==0:
            result = "Completed"
        elif error>0 and error<len(ERRORS):
            result = "ERROR: "+ERRORS[job.error]
        else:
            result = "ERROR: "+str(job.error)
        #assert False
    else:
        result = "In progress"
    return result

    