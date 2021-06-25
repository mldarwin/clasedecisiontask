# DOSE Code

Code to implement the Dynamically Optimized Sequential Experimentation (DOSE) approach to optimizing the order of choice presentation in a risky decision-making framework such as CLASE.

See _Wang, Filiba, & Camerer (2010)_, and _Chapman, Snowberg, Wang, & Camerer (2019)_ -> [link to the paper](http://jnchapman.com/assets/pdf/dose.pdf).

From June 2021 email communication with Dr. Wang:

> Please see some code attached - the procedure calls various functions and these are in separate files. The file that runs everything is called runDose.m, then the other files call the component parts of the procedure. Let me know if you run into any problems.
>
> You can use the "runDose" function  to run the procedure directly from matlab - you just need to specify an input file which has choice data (I'm including an example of the format used) and then information about the priors you want to use. The information at the top of the .m should tell you most of what you need to know. For 'utilityFunction' you should enter 'CRRA' for it to work with this code. We have found a uniform discretized prior works well - you can enter the number of points to use in the runDose function. Alternatively, you can load a custom file with the details of the specific prior you need.
>
> Here is an example of the command I am using, in case it helps: runDose('EconoSurvRiskPrefWave1','CRRA',100,100,100,'original','WEPwave1_CRRAoriginal100pt','maximum','allSubjects','no', 'header')
>
> If you just want to get a sense of how we got this to work, the most relevant files would be shChoiceOptimal and KLCalculator.

From June 2021 email communication with Dr. Chapman:

> That code is organized to identify a question sequence where we already have a set of answers. i.e., it does not gather information dynamically, although it should be straightforward to reorganize it to do so. However, if you are only planning a short question sequence, it may be best to implement the procedure the same way we did in the survey: identify all possible response sets and corresponding question sequences in advance. That could be helpful in that it avoids processing time between each question round (while the procedure identifies the next optimal question).

... and...

> We have created two different question sequences: one was 10Q (so around 1000 possible routes) and the latter 20Q (so around 1million). The former ("simple" model) simply allowed for risk aversion and loss aversion parameters (and an error parameter), the latter ("complex" model) allowed for separate risk aversion parameters over gains or losses.
>
> Using a reasonably high-powered machine (16GB Ram, 4 cores) the 10Q sequence was very quick, the 20Q was more burdensome. From memory it took 2-3 days, splitting the estimation into multiple parts so I could run them simultaneously.  
>
> With the simple model, my sense is that more than 20Q would be possible, but with the complex model it would be difficult at least without turning to a server allowing a very large amount of parallel processing.
>
> Re the time for optimal question selection: the new version of the code is very quick - probably less than a second for the simple model, a few seconds for the complex one. (Of course, it would take more time to allow for the data to be inputted dynamically as well). The time taken really varies according to the size of the question set you are using and how refined a prior you use for the question selection procedure. So I do think dynamic question selection is very feasible - my only concern would be how it performs on a lower spec machine (I simply haven't tried).
>
> For the 10Q procedure, which is reported in the current version of our paper (here)[http://jnchapman.com/assets/pdf/dose.pdf], the question set was around 1,700 questions - you can see some details in online appendix B.3). For the newer module we had to reduce the question set to accommodate the greater dimensionality of the prior (because of the additional parameter).
