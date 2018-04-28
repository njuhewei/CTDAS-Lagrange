#! /usr/bin/env python
# rc.py


# ------------------------------------------------
# help
# ------------------------------------------------

"""
Deal with model settings in `rc` format.

RCFILES

    A rcfile is a text file with key/value pairs seperated by a ':', e.g.

      my.flag    :  T
      my.answer  :  42

    The following functionality is supported:

     * Empty lines are ignored.

     * Comment lines are introduced by a '!' as first character.

     * Long values could be continued at the next line after a '\' as last character.

     * Include the key/value pairs from another file:

         #include an/other.rc

     * Substitute environment variables when available:

         tm.dir : ${HOME}/TM5/cy3

     * Substitute values of other keys in a line:

         build.dir            :  ${tm.dir}/build
         grid                 :  glb300x200
         input.${grid}.path   :  /data/input/${grid}

       Substitions are allowed in both key names as well as values.
       The substitutions are performed in a loop until nothing
       has to be substituted anymore, or some substitutions could
       not be applied at al; for the later an error is raised.
       Values to be substituted could therefore be set before and
       after they are used.

       Note that if a key has the same name as an environment variable,
       the new value will be assigned to the key instead of the value
       retrieved from the environment:

         HOME      :  /some/other/dir/

     * Substitude some specials:

         ${pid}     # evaluates to the current process id; 
                    # useful for names of log files etc
         ${script}  # evaluates to the base name of the calling script, 
                    # thus without .py etc
                    
     * Instead of variables of the form '${..}' other patterns could be
       specified with the optional 'marks' tupple (see below).

     * Old-style '#eval' lines are still supported:

         #eval RUNDIR = /path/to/mydir
         tmdir : ${RUNDIR}/TM5/cy3

       In this example, the value of RUNDIR will be evaluated and substituted 
       in all {key,value} pairs. This feature is obsolete and a warning will 
       be issued. The proper way to use this is with {key,value} pairs too:

         run.dir   : /path/to/mydir
         tmdir     : ${run.dir}/TM5/cy3
         
     * Comment starting with '!' is stripped from the values.
       To have a value including exclamation marks, use '\!' but do
       not expect that the rest of the value is scanned for comment too:
       
           my.value      :   -999    ! just an integer value
           my.message    :   This value has 64 characters \! Count if you don't believe it ...

     * If you trust yourself you might try to use conditional expressions:
     
           #if ${my.number} == 1
           message    : Welcome
           #else
           message    : Whatever ...
           #endif
           
       The conditions should be valid python expressions that evaluate to a boolean;
       value substitutions are performed before evaluation. Examples:

                ${my.runmode} == 4
                "${my.tracer}" == "CH4"

       Keep it simple! Very complicated and nested if-statements might not be
       resolved correctly, and are in any case not easy to understand for other users!
       
       In the example above, an exception could be raised by the special error expression;
       everything behind the '#error' mark is displayed as an error message:
       
            #error No settings provided for value : ${my.value}


USAGE AS SCRIPT

    Called in script form, the following syntaxis is supported:
    
       rc.py [options] <rcfile> <key>
       rc.py -h|--help
       
    The <rcfile> is read and the value defined by <key> is printed
    to the standard output.

    Use the --help option for more documentation.

    
USAGE AS PYTHON MODULE

    Import the module with:
    
        import rc

    Initialiase by reading all settings in a rcfile,
    supporting the functionality described in the 'RCFILES' section.

        rcf = RcFile( 'settings.rc' )

    The initialisation accepts some optional arguments.
    Set the silent flag to True to ignore warnings.

        rcf = RcFile( 'settings.rc', silent=False )

    Use the optional 'marks' tupple to define that variables to be expanded
    are marked other than '${..}' but rather '<mark1>..<mark2>' :

        rcf = RcFile( 'settings.rc', marks=('${','}') )

    Test to see if a key is defined:

        if rcf.has_key('my.flag') :
            print 'value of my flag is : ', rcf['my.flag']
    
    Extract a list with all keys:
    
        rcf.keys()

    A 'get' function is provided to extract values:

     * by default, the 'get' function returns the value as a str type:

         s = rcf.get('my.value')
    
     * a second argument is the name of the python type to which
       the value is converted to:
    
         i = rcf.get('my.flag','int')

     * if the return value should be a 'bool', the result is
         True  for values     : 'True' , 'T', 'yes', or '1' ,
         and False for value  : 'False', 'F', 'no' , or '0' ;
       for other values an error is raised;
     
     * return a default value if the key is not found:

            rcf.get( 'my.flag', default=False )
    
     * print a debug message to the logging system for each extracted key:
     
            rcf.get( 'my.flag', verbose=True ) 

    Add a new value, comment is optional:

        rcf.add( 'my.iter', 2, comment='iteration number for restart' )

    Assign a new value to an existing key:

        rcf.replace( 'my.flag', True )

    Scan a character line for all occurances of ${<key>} and subsitute for
    the rc value assigned to <key> :

        line = rcf.substitute( line )

    Write the dictionary (with all variables expanded and included files included)
    to new file:

         rcf.write('newfile.rc')

         
USAGE AS PYTHON MODULE - BACKWARDS COMPATIBILITY

    For backwards compatibility with older implementations of the rc.py module,
    two extra routines are available.
    
    To read rc-file by making an instance of the RcFile class, 
    and to returns a dictionary of values only, use:
            
        rcdict = read( 'test.rc' [,silent=False] )
        
    Create a new rcfile and fill with key/values from a dictionary:

        write( 'test.rc', rcdict )
        

HISTORY

    2008? Andy Jacobson, NOAA
      Translation to python of original shell script 'go_readrc' .
    2009-06 Wouter Peters, WUR
      Support substitution of previously defined variables.
    2009-06 Arjo Segers, TNO
      Support include files.
    2009-09 Arjo Segers, TNO
      Re-coded into class.
      Implemented substitution loop.
    2009-11 Arjo Segers, JRC
      Added main program to run this file as a shell script.
      Added replace and substitute routines.
    2010-03 Arjo Segers, JRC
      Support simple if-statements.
      Support comment in values.
    2010-07 Wouter Peters, WUR
      Downgraded to work for python 2.4.3 too.
      Added read/write routines for backwards compatibility.
    2010-07-27 Arjo Segers, JRC
      Maintain list with rcfile names and line numbers to be displayed
      with error messages to identify where problematic lines are found.
    2010-07-28 Andy Jacobson, NOAA
      Add second dictionary of key,linetrace values to help track the 
      provenance of #included keys (to debug multiple key instances).
      Identify duplicate keys by checking on different source lines
      instead of checking if the values are different.
"""
import re
import os
import sys
import logging
        

# ------------------------------------------------
# classes
# ------------------------------------------------


class RcFile(object) :

    """
    Class to store settings read from a rcfile.
    """   

    def __init__(self, filename, silent=False, marks=('${', '}')) :

        """ 
        
        Usage:
        
          rcf = RcFile( 'settings.rc' [,silent=False] [marks=('${','}')] )

        Read an rc-file and expand all the keys and values into a dictionary.
        Do not shout messages if silent is set to True. 
        The 2-item tupple (mark1,mark2) could be used to re-define the default
        key pattern '${..}' into something else:
          <mark1>...<mark2>

        """

        # info ...
        logging.debug('reading rcfile %s ...' % filename)

        # check ...
        if not os.path.exists(filename) :
            msg = 'rcfile not found : %s' % filename ; logging.error(msg)
            raise IOError, msg
        #endif
        
        # store file name:
        self.filename = filename
        # store rc-file root directory:
        self.rootdir = os.path.split(filename)[0]
        
        # storage for processed rcfile:
        self.outfile = []
        
        # storage for key/value pairs:
        self.values = {}

        # storage for key/source file pairs:
        self.sources = {}
        
        # open the specified rc-file:
        f = open(filename, 'r')
        # store all lines in a list:
        inpfile = f.readlines()
        # close:
        f.close()
        
        # create traceback info:
        inptrace = []
        for jline in range(len(inpfile)) :
            inptrace.append('"%s", line %-10s' % (filename, str(jline + 1)))
        #endfor

        # flags:
        warned_for_eval = False
        
        # pass counter:
        ipass = 1
        
        # loop until all substitutions and inclusions are done:
        while True :
        
            # start again with empty output file:
            self.outfile = []
            # also empty traceback info:
            self.trace = []
            # init current line:
            line = ''
            # assume nothing has to be done after this loop:
            something_done = False
            something_to_be_done = False
            # maintain list with unresolved lines (with keys that could not be evaluated yet):
            unresolved_lines = []
            # maintain list with keys from which the value could not be resolved yet:
            keys_with_unresolved_value = []
            # maintain list with undefined keys; 
            # some might be part of the keys_with_unresolved_value list:
            undefined_keys = []
                        
            # stack for conditional evaluation;
            # each element is a tuple with elements:
            #   resolved (boolean) true if the if-statement is evaluated
            #   flag     (boolean) true if the lines below the statement 
            #              are to be included
            #   anyflag    (boolean) used to check if any of the 'if' or 'elif' conditions
            #              in this sequence evaluated to True
            #   line     (char) description of the line for messages
            ifstack = []

            #print ''
            #print '---[pass %i]-------------------------------------' % ipass
            #for line in inpfile : print line.strip()
            
            # loop over lines in input file:
            iline = -1
            for inpline in inpfile :
            
                # line counter:
                iline = iline + 1
                
                # cut current traceback info from list:
                linetrace_curr = inptrace.pop(0)
                
                # set full traceback info:
                if line.endswith('\\') :
                    # current input line is a continuation; combine:
                    qfile, qlinenrs = linetrace.split(',')
                    qnrs = qlinenrs.replace('lines', '').replace('line', '')
                    if '-' in qnrs :
                        qnr1, qnr2 = qnrs.split('-')
                    else :
                        qnr1, qnr2 = qnrs, qnrs
                    #endif
                    linetrace = '%s, lines %-9s' % (qfile, '%i-%i' % (int(qnr1), int(qnr2) + 1))
                else :
                    # just copy:
                    linetrace = linetrace_curr
                #endif
                    
                # remove end-of-line character:
                inpline = inpline.strip()
                
                ## DEBUG: display current line ...
                #print '%4i | %s' % (iline,inpline)
                #print '%4i | %s     %s' % (iline,inpline,linetrace)
                
                #
                # empty lines
                #

                # skip empty lines:
                if len(inpline) == 0 :
                    # add empty line to output:
                    self.outfile.append('\n')
                    # add traceback info:
                    self.trace.append(linetrace)
                    # next will be a new line:
                    line = ''
                    # next input line:
                    continue
                #endif
                
                #
                # comment lines
                #

                # skip comment:
                if inpline.startswith('!') :
                    # add copy to output file:
                    self.outfile.append('%s\n' % inpline)
                    # add traceback info:
                    self.trace.append(linetrace)
                    # next will be a new line:
                    line = ''
                    # next input line:
                    continue
                #endif
                
                #
                # continuation lines
                #

                # current line has continuation mark '\' at the end ?
                # then add this input line:
                if line.endswith('\\') :
                    # remove continuation character, add input line:
                    line = line[:-1] + inpline
                else :
                    # bright new line:
                    line = inpline
                #endif

                # continuation mark ? then next line of input file:
                if line.endswith('\\') : continue
                
                #
                # conditional settings (1)
                #

                ## debug ...
                #print 'xxx0 ', ifstack
                
                # is this the begin of a new condition ?
                mark = '#if'
                if line.startswith(mark) :
                    # push temporary flag to stack, will be replaced after evaluation of condition:
                    ifstack.append((False, True, False, linetrace))
                    # debug ...
                    #print 'xxx1 ', ifstack
                #endif

                mark = '#elif'
                if line.startswith(mark) :
                    # check ...
                    if len(ifstack) == 0 :
                        logging.error('found orphin "%s" in %s' % (mark, linetrace))
                        raise Exception
                    #endif
                    # remove current top from stack:
                    resolved, flag, anyflag, msg = ifstack.pop()
                    # did one of the previous #if/#elif evaluate to True already ?
                    if resolved and anyflag :
                        # push to stack that this line resolved to False :
                        ifstack.append((True, False, anyflag, linetrace))
                        # copy to output:
                        self.outfile.append('%s\n' % line)
                        # add traceback info:
                        self.trace.append(linetrace)
                        # next input line:
                        continue
                    else :
                        # push temporary flag to stack, will be replaced after evaluation of condition:
                        ifstack.append((False, True, anyflag, linetrace))
                    #endif
                    ## debug ...
                    #print 'xxx2 ', ifstack
                #endif

                mark = '#else'
                if line.startswith(mark) :
                    # check ...
                    if len(ifstack) == 0 :
                        logging.error('found orphin "%s" in %s' % (mark, linetrace))
                        raise Exception
                    #endif
                    # remove current top from stack:
                    resolved, flag, anyflag, msg = ifstack.pop()
                    # get higher level settings:
                    if len(ifstack) > 0 :
                        resolved_prev, flag_prev, anyflag_prev, msg_prev = ifstack[-1]
                    else :
                        flag_prev = True
                    #endif
                    # should next lines be included ? 
                    # reverse flag, take into acount higher level:
                    new_flag = (not flag) and (not anyflag) and flag_prev
                    # push to stack:
                    ifstack.append((resolved, new_flag, False, linetrace))
                    # debug ...
                    #print 'xxx3 ', ifstack
                    # copy to output:
                    self.outfile.append('%s\n' % line)
                    # add traceback info:
                    self.trace.append(linetrace)
                    # next input line:
                    continue
                #endif

                # is this the end of a condition ?
                mark = '#endif'
                if line.startswith(mark) :
                    # check ...
                    if len(ifstack) == 0 :
                        logging.error('found orphin "%s" in %s' % (mark, linetrace))
                        raise Exception
                    #endif
                    # remove top from stack:
                    top = ifstack.pop()
                    # copy to output:
                    self.outfile.append('%s\n' % line)
                    # add traceback info:
                    self.trace.append(linetrace)
                    # next input line:
                    continue
                #endif

                # within if-statements ?
                if len(ifstack) > 0 :
                    # extract current top:
                    resolved, flag, anyflag, msg = ifstack[-1]
                    # already resolved ? then check if this line should be skipped:
                    if resolved and (not flag) :
                        # skip this line, but include commented version in output:
                        self.outfile.append('!%s\n' % line)
                        # add traceback info:
                        self.trace.append(linetrace)
                        # read the next input line:
                        continue
                    #endif
                #endif
                
                #
                # handle '#eval' lines
                #

                mark = '#eval'        
                if line.startswith(mark):
                    # info ..
                    if not warned_for_eval :
                        if not silent: logging.warning('the #eval statements in rc-files are deprecated, use {key:value} pairs instead')
                        warned_for_eval = True
                    #endif
                    # add commented copy to output:
                    self.outfile.append('!evaluated>>> ' + line)
                    # add traceback info:
                    self.trace.append(linetrace)
                    # remove leading mark:
                    line = line.lstrip(mark)
                    # multiple settings are seperated by ';' :
                    evals = line.split(';')
                    # insert in output file:
                    for k in range(len(evals)) :
                        # split in key and value:
                        key, value = evals[k].split('=')
                        # insert:
                        self.outfile.append('%s : %s' % (key, value))
                        # add traceback info:
                        self.trace.append(linetrace)
                    #endfor
                    # next input line:
                    continue
                #endif

                #
                # replace ${..} values
                #
                
                # ensure that common marks are evaluated correctly:
                start_mark = marks[0].replace('{', '\{').replace('<', '\<').replace('$', '\$')
                close_mark = marks[1].replace('}', '\}').replace('>', '\>')
        
                # set syntax of keywords to be matched, e.g. '${...}' :
                pattern = start_mark + '[A-Za-z0-9_.]+' + close_mark

                # make a regular expression that matches all variables:
                rc_varpat = re.compile(pattern)

                # search all matching paterns:
                pats = re.findall(rc_varpat, line)
                # counter for unexpanded substitutions:
                ntobedone = 0
                # loop over matches:
                for pat in pats :
                    # remove enclosing characters:
                    key = pat.lstrip(start_mark).rstrip(close_mark)
                    # test some dictionaries for matching key:
                    if self.values.has_key(key) :
                        # get previously defined value:
                        val = self.values[key]
                        # substitute value:
                        line = line.replace(pat, val)
                        # set flag:
                        something_done = True
                    elif os.environ.has_key(key) :
                        # get value from environment:
                        val = os.environ[key]
                        # substitute value:
                        line = line.replace(pat, val)
                        # set flag:
                        something_done = True
                    elif key == 'pid' :
                        # special value: process id; convert to character:
                        val = '%i' % os.getpid()
                        # substitute value:
                        line = line.replace(pat, val)
                        # set flag:
                        something_done = True
                    elif key == 'script' :
                        # special value: base name of the calling script, without extension:
                        script, ext = os.path.splitext(os.path.basename(sys.argv[0]))
                        # substitute value:
                        line = line.replace(pat, script)
                        # set flag:
                        something_done = True
                    else :
                        # could not substitute yet; set flag:
                        ntobedone = ntobedone + 1
                        # add to list with unresolved keys:
                        if key not in undefined_keys : undefined_keys.append(key)
                        # continue with next substitution:
                        continue
                    #endif
                #endfor  # matched patterns
                # not all substituted ?
                if ntobedone > 0 :
                    # add line to output:
                    self.outfile.append(line)
                    # add traceback info:
                    self.trace.append(linetrace)
                    # a new pass is needed:
                    something_to_be_done = True
                    # store for info messages:
                    unresolved_lines.append('%s | %s' % (linetrace, line))
                    # could this be a 'key : value' line ?
                    if ':' in line :
                        # split, remove leading and end space:
                        qkey, qvalue = line.split(':', 1)
                        qkey = qkey.strip()
                        # assume it is indeed a key if it does not contain whitespace,
                        # no start mark, and does not start wiht '#' :
                        if (' ' not in qkey) and (start_mark not in qkey) and (not qkey.startswith('#')) :
                            # add to list:
                            if qkey not in keys_with_unresolved_value : keys_with_unresolved_value.append(qkey)
                        #endif
                    # next input line:
                    continue
                #endif

                #
                # handle '#include' lines
                #

                mark = '#include'
                if line.startswith(mark) :
                    # remove leading mark, what remains is file to be included:
                    inc_file = line.lstrip(mark).strip()
                    # check ...
                    if not os.path.exists(inc_file) :
                        inc_file = os.path.join(self.rootdir, inc_file)
                        logging.debug('Added rootdir to requested include: %s' % inc_file)
                    #endif
                    if not os.path.exists(inc_file) :
                        logging.error('include file not found : %s' % inc_file)
                        logging.error(linetrace)
                        raise IOError, 'include file not found : %s' % inc_file
                    #endif
                    # read content:
                    inc_f = open(inc_file, 'r')
                    inc_rc = inc_f.readlines()
                    inc_f.close()
                    # add extra comment for output file:
                    self.outfile.append('! >>> %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n' % inc_file)
                    self.outfile.extend(inc_rc)
                    self.outfile.append('! <<< %s <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n' % inc_file)
                    # add traceback info:
                    self.trace.append(linetrace)
                    for jline in range(len(inc_rc)) :
                        self.trace.append('"%s", line %-10s' % (inc_file, str(jline + 1)))
                    #endfor
                    self.trace.append(linetrace)
                    # set flag:
                    something_done = True
                    # a new pass is needed:
                    something_to_be_done = True
                    # next input line:
                    continue
                #endif


                #
                # conditional settings (2)
                #

                # evaluate conditional expressions:
                mark1 = '#if'
                mark2 = '#elif'
                if line.startswith(mark1) or line.startswith(mark2) :
                    # remove leading mark, what remains is logical expression:
                    expression = line.lstrip(mark1).strip()
                    expression = line.lstrip(mark2).strip()
                    # common mistake is to add a ':' as in python; remove this:
                    if expression.endswith(':') : expression = expression.rstrip(':').strip()
                    # evaluate:
                    try :
                        flag = eval(expression)
                    except :
                        logging.error('could not evaluate expression:')
                        logging.error('    %s' % expression)
                        logging.error('in %s' % linetrace)
                        raise Exception
                    #endtry
                    # remove temporary top added before during this pass:
                    tmp_statement, tmp_flag, tmp_anyflag, tmp_msg = ifstack.pop()
                    # extract current top if necessary:
                    if len(ifstack) > 0 :
                        dummy_statement, prev_flag, dummy_anyflag, dummy_msg = ifstack[-1]
                    else :
                        prev_flag = True
                    #endif
                    # should next lines be included ?
                    new_flag = prev_flag and tmp_flag and flag
                    # any if/elif evaluated to true in this sequence ?
                    new_anyflag = tmp_anyflag or new_flag
                    # add to stack, now resolved, take into accout current flag:
                    ifstack.append((True, new_flag, new_anyflag, linetrace))
                    # debug ...
                    #print 'xxx2 ', ifstack
                    # copy to output:
                    self.outfile.append('%s\n' % line)
                    # add traceback info:
                    self.trace.append(linetrace)
                    # next input line:
                    continue
                #endif

                #
                # error message
                #
                
                # special command to rais an exception:
                mark = '#error'
                if line.startswith(mark) :
                    # remove leading mark, what remains is error message:
                    msg = line.lstrip(mark).strip()
                    # display:
                    logging.error(msg)
                    # add info:
                    logging.error('error message in %s' % linetrace)
                    # stop:
                    raise Exception
                #endif

                #
                # checks
                #

                # common mistake ...
                if line.startswith('#') :
                    logging.error('line in rcfile starts with "#" but is not an "#include" or other special line;')
                    logging.error('if it is supposed to be comment, please start with "!" ...')
                    logging.error('  %s' % line)
                    logging.error('%s' % linetrace)
                    raise IOError
                #endif

                # check ...
                if ':' not in line :
                    logging.error('key/value line should contain a ":"')
                    logging.error('%s' % linetrace)
                    raise IOError
                #endif

                #
                # add to output
                #

                # add line to output:
                self.outfile.append('%s\n' % line)
                # add traceback info:
                self.trace.append(linetrace)

                #
                # add key/value pair
                #
                
                # not if inside an unresolved if-statement ...
                if len(ifstack) > 0 :
                    # get top values:
                    resolved, flag, anyflag, msg = ifstack[-1]
                    # not resolved yet ? then continue:
                    if not resolved : continue
                #endif
                
                # split in key and value; 
                # value might contain ':' too, so at maximum 1 split:
                key, val = line.split(':', 1)
                
                # remove comment from value:
                if '!' in val :
                    # not if '\!' is in the value ...
                    if not '\!' in val : val, comment = val.split('!')
                    # replace all slash-comments:
                    val = val.replace('\!', '!')
                #endif

                # remove spaces:
                key = key.strip()
                val = val.strip()

                # already defined ?
                if self.values.has_key(key) :
                    # this will occure often after the first pass since
                    # the keys are resolved again and again ;
                    # therefore, only complain if this definition is read
                    # from a different line :
                    if linetrace != self.sources[key] :
                        logging.error('duplicated key \"%s\" found:' % key)
                        logging.error('first definition in %s is:' % self.sources[key])
                        logging.error('  %s  : %s' % (key, str(self.values[key])))
                        logging.error('second definition in %s is:' % linetrace.strip())
                        logging.error('  %s  : %s' % (key, str(val)))
                        raise Exception
                    #endif
                else :
                    # store new value:
                    self.values[key] = val
                    self.sources[key] = linetrace
                    # set flag:
                    something_done = True
                #endif

                # display key and value ...
                #print '                                --> %s : %s, from %s' % (key,val, linetrace)

            #endfor  # loop over lines in text
            
            ## info ...
            #print '~~~ outfile ~~~~~~~~~~~~~~~~~~~~~~~'
            #for line in self.outfile : print line.strip()
            #print '~~~ key/values ~~~~~~~~~~~~~~~~~~~~'
            #for k,v in self.iteritems() :
            #    print '%s  :  %s' % (k,v)
            ##endfor
            #print '-------------------------------------------------'
            #print ''
            
            # check ...
            if len(ifstack) > 0 :
                logging.error('unterminated if-statement ; current stack:')
                for resolved, flag, anyflag, msg in ifstack : logging.error(msg)
                logging.error('please fix the rcfile or debug this script ...')
                raise Exception
            #endif

            # check ...
            if something_to_be_done :
                # check for unterminated loop ...
                if not something_done :
                    # list all unresolved lines:
                    logging.error('Could not resolve the following lines in rcfile(s):')
                    logging.error('')
                    for uline in unresolved_lines :
                        logging.error('    %s' % uline)
                    #endfor
                    logging.error('')
                    # list all undefined keys:
                    logging.error('  Undefined key(s):')
                    logging.error('')
                    for ukey in undefined_keys :
                        # do not list them if they are undefined because the value
                        # depends on other undefined keys:
                        if ukey not in keys_with_unresolved_value :
                            # display:
                            logging.error('    %s' % ukey)
                            # loop over unresolved lines to see in which the key is used:
                            for uline in unresolved_lines :
                                # search for  '${key}' pattern:
                                if marks[0] + ukey + marks[1] in uline :
                                    logging.error('      %s' % uline)
                                #endif
                            #endfor
                            logging.error('')
                        #endif
                    #endfor
                    logging.error('please fix the rcfile(s) or debug this script ...')
                    raise Exception
                #endif
            else :
                # finished ...
                break
            #endif
            
            # for safety ...
            if ipass == 100 :
                logging.error('resolving rc file has reached pass %i ; something wrong ?' % ipass)
                raise Exception
            #endif
            
            # new pass:
            ipass = ipass + 1
            # renew input:
            inpfile = self.outfile
            # renew traceback:
            inptrace = self.trace
            
        #endwhile   # something to be done
        
    #enddef  # __init__
    
    
    # ***
    
    
    def has_key(self, key) :
    
        # from dictionairy:
        return self.values.has_key(key)
        
    #enddef
    
    
    # ***
    
    
    def keys(self) :
    
        # from dictionairy:
        return self.values.keys()
        
    #enddef
    
    
    # ***


    def get(self, key, totype='', default=None, verbose=False) :
    
        """
        rcf.get( 'my.value' [,default=None] )
        Return element 'key' from the dictionairy.
        If the element is not present but a default is specified, than return
        the default value.
        If 'verbose' is set to True, then print debug messages to the logging
        about which values is returned for the given key.
        The option argument 'totype' defines the conversion to a Python type.
        If 'totype' is set to 'bool', the return value is the
        boolean True for values 'T', 'True', 'yes', and '1',
        and False for 'F', 'False', 'no', or '0' ;
        for other values, an exception will be raised.
        """
        
        # element found ?
        if self.values.has_key(key) :
            # copy value:
            value = self.values[key]
            # convert ?
            if totype == 'bool' :
                # convert to boolean:
                if value in ['T', 'True', 'yes', '1'] :
                    value = True
                elif value in ['F', 'False', 'no', '0'] :
                    value = False
                else :
                    logging.error("value of key '%s' is not a boolean : %s" % (key, str(value)))
                    raise Exception
                #endif
            elif len(totype) > 0 :
                # convert to other type ...
                value = eval('%s(%s)' % (totype, value))
            #endif
            # for debugging ...
            if verbose : logging.debug('rc setting "%s" : "%s"' % (key, str(value)))
        else :
            # default value specified ?
            if default != None :
                # copy default:
                value = default
                # for debugging ...
                if verbose : logging.debug('rc setting "%s" : "%s" (deault)' % (key, str(value)))
            else :
                # something wrong ...
                logging.error("key '%s' not found in '%s' and no default specified" % (key, self.filename))
                raise Exception
            #endif
        #endif
        
        # ok
        return value
        
    #enddef
    
    
    # ***
    
    
    def replace(self, key, val) :
    
        """
        Replace a key by a new value.
        """
        
        # search for a line '<key>   : <val>' 
        # loop over lines in output file:
        found = False
        for iline in range(len(self.outfile)) :
            # extract:
            line = self.outfile[iline]
            # skip lines that are no key:value pair for sure ...
            if ':' not in line : continue
            # split once at first ':'
            k, v = line.split(':', 1)
            # match ?
            if k.strip() == key :
                # replace line in original file:
                self.outfile[iline] = '%s : %s\n' % (k, str(val))
                # replace value:
                self.values[key] = val
                # set flag:
                found = True
                # found, thus no need to continue:
                break
            #endif
        #endfor  # lines
        # not found ?
        if not found :
            logging.error('could not replace key : %s' % key)
            raise Exception
        #endif
        
        # ok
        return
    
    #enddef
    
    
    # ***
    
    
    def add(self, key, val, comment='') :
    
        """Add a new key/value pair."""
        
        # add lines:
        self.outfile.append('\n')
        if len(comment) > 0 : self.outfile.append('! %s\n' % comment)
        self.outfile.append('%s : %s\n' % (key, str(val)))

        # add to dictionairy:
        self.values[key] = val
        
        # ok
        return
    
    #enddef
    
    
    # ***
    
    
    def substitute(self, line, marks=('${', '}')) :
    
        """
        Return a line with all '${..}' parts replaced by the corresponding rcfile values.
        The 2-item tupple (mark1,mark2) could be used to re-define the default
        key pattern '${..}' into something else:
          <mark1>...<mark2>
        """
        
        # ensure that common marks are evaluated correctly:
        start_mark = marks[0].replace('{', '\{').replace('<', '\<').replace('$', '\$')
        close_mark = marks[1].replace('}', '\}').replace('>', '\>')

        # set syntax of keywords to be matched, e.g. '${...}' :
        pattern = start_mark + '[A-Za-z0-9_.]+' + close_mark

        # make a regular expression that matches all variables:
        rc_varpat = re.compile(pattern)

        # search all matching paterns:
        pats = re.findall(rc_varpat, line)
        # loop over matches:
        for pat in pats :
            # remove enclosing characters:
            key = pat.lstrip(start_mark).rstrip(close_mark)
            # test dictionary for matching key:
            if self.values.has_key(key) :
                # get previously defined value:
                val = self.values[key]
                # substitute value:
                line = line.replace(pat, val)
            #endif
        #endfor  # matched patterns

        # ok
        return line
        
    #enddef
    

    # ***


    def WriteFile(self, filename) :

        """ write the dictionary to file"""

        # open file for writing:
        f = open(filename, 'w')

        ## loop over key/value pairs:
        #for k,v in self.iteritems():
        #    # add line; at least the specified number of characters 
        #    # is used for the key:
        #    f.write( '%-20s:%s\n' % (k,v) )
        ##endfor

        # write processed input:
        f.writelines(self.outfile)
        
        # close file:
        f.close()
        
    #endif
    

#endclass    # RcFile


# ***


def read(rcfilename, silent=False) :

    """ 
    This method reads an rc-file by making an instance of the RcFile class, 
    and then returns the dictionary of values only. 
    This makes it backwards compatible with older implementations of the rc.py module
    """

    rcdict = RcFile(rcfilename, silent=silent)

    return rcdict.values

#enddef


# ***


def write(filename, rcdict) :

    """
    This method writes an rc-file dictionary. 
    This makes it backwards compatible with older implementations of the rc.py module
    """

    # open file for writing:
    f = open(filename, 'w')

    # loop over key/value pairs:
    for k, v in rcdict.items():
        # add line; at least the specified number of characters 
        # is used for the key:
        f.write('%-20s:%s\n' % (k, v))
    #endfor

    # close file:
    f.close()

#enddef



# ------------------------------------------------
# script
# ------------------------------------------------


if __name__ == '__main__':

    # external ...
    import optparse
    import traceback
    
    # extract arguments from sys.argv array:
    #   0 = name of calling script, 1: = actual arguments
    args = sys.argv[1:]
    
    # set text for 'usage' help line:
    usage = "\n    %prog <rcfile> <key> [-b|--bool] [--default<=value>]\n    %prog <rcfile> -w|--write\n    %prog -h|--help\n    %prog -d|--doc"

    # initialise the option parser:
    parser = optparse.OptionParser(usage=usage)
    
    # define options:
    parser.add_option("-d", "--doc",
                         help="print documentation",
                         dest="doc", action="store_true", default=False)
    parser.add_option("-v", "--verbose",
                         help="print information messages",
                         dest="verbose", action="store_true", default=False)
    parser.add_option("-b", "--bool",
                         help="return 'True' for values 'T', 'True', 'yes', or '1', and 'False' for 'F', 'False', 'no', or '0'",
                         dest="boolean", action="store_true", default=False)
    parser.add_option("--default",
                         help="default value returned if key is not found",
                         dest="default", action="store")
    parser.add_option("-w", "--write",
                         help="write pre-processed rcfile",
                         dest="write", action="store_true", default=False)
    
    # now parse the actual arguments:
    opts, args = parser.parse_args(args=args)
    
    # print documentation ?
    if opts.doc :
        print __doc__
        sys.exit(0)
    #endif
    
    # recfile argument should be provided:
    if len(args) < 1 :
        parser.error("no name of rcfile provided\n")
    #endif
    # extract:
    rcfile = args[0]
    
    # read rcfile in dictionary:
    try :
        rcf = RcFile(rcfile, silent=(not opts.verbose))
    except :
        if opts.verbose : logging.error(traceback.format_exc())
        sys.exit(1)
    #endtry
    
    # print pre-processed file ?
    if opts.write :
        for line in rcf.outfile : print line.strip()
        sys.exit(0)
    #endif

    # key argument should be provided:
    if len(args) < 2 :
        parser.error("no name of rckey provided\n")
    #endif
    # extract:
    rckey = args[1]
    
    # key present ?
    if rcf.has_key(rckey) :

        # print requested value:
        if opts.boolean :
            # extract value:
            flag = rcf.get(rckey, 'bool')
            # print result:
            if flag :
                print 'True'
            else :
                print 'False'
            #endif
        else :
            # extract value:
            value = rcf.get(rckey)
            # display:
            print value
        #endif
        
    else :

        # default value provided ?
        if opts.default != None :
            # display:
            print opts.default
        else :
            print 'ERROR - key "%s" not found in rcfile "%s" and no default specified' % (rckey, rcfile)
            sys.exit(1)
        #endif

    #endif
    
#endif


# ------------------------------------------------
# end
# ------------------------------------------------
