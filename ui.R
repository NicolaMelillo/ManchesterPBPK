shinyUI(navbarPage("Manchester PBPK - beta",
                   tabPanel("Inputs/Simulation",
                            sidebarLayout(
                              sidebarPanel(
                                helpText("Input key paramaters and press run simulation."),
                                fluidRow(
                                  column(6,
                                         h5('Route')
                                  ),
                                  column(6,
                                         radioButtons("Route", 
                                                      label = NULL, 
                                                      choices = list("Bolus iv" = 1, 
                                                                     "Bolus po solid" = 2,
                                                                     "Bolus po dissolved" = 3),
                                                      selected = 1)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Schedule')
                                  ),
                                  column(6,
                                         radioButtons("Schedule", 
                                                            label = NULL, 
                                                            choices = list("Once-Daily" = 1, 
                                                                           "Twice-Daily" = 2),
                                                            selected = 1)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Dose (mg)')
                                  ),
                                  column(6,
                                         numericInput("dose", label = NULL, value = 10)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Days')
                                  ),
                                  column(6,
                                         numericInput("days", label = NULL, value = 1, min=0, max=30, step=1)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('MW (g/Mol)')
                                  ),
                                  column(6,
                                         numericInput('MW', label = NULL, value = 500)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Type')
                                  ),
                                  column(6,
                                         radioButtons("type", 
                                                            label = NULL, 
                                                            choices = list("Neutral" = 0, 
                                                                           "Acid" = 1,
                                                                           "Base" = 2),
                                                            selected = 0)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('logPow')
                                  ),
                                  column(6,
                                         numericInput('logPow', label = NULL, value = 2)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Fup')
                                  ),
                                  column(6,
                                         numericInput('fup', label = NULL, value = 0.8,min=0,max=1)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('BP')
                                  ),
                                  column(6,
                                         numericInput('BP', label = NULL, value = 0.8)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('pKa')
                                  ),
                                  column(6,
                                         numericInput('pKa', label = NULL, value = 0.8)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Particle Size (um)')
                                  ),
                                  column(6,
                                         numericInput('r', label = NULL, value = 2.5)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Density of Formulation (g/L)')
                                  ),
                                  column(6,
                                         numericInput('rho', label = NULL, value = 1000)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Intrinisic Solubility (mg/L)')
                                  ),
                                  column(6,
                                         numericInput('Csint', label = NULL, value = 100)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Effective Permeability (10^-4 cm/s)')
                                  ),
                                  column(6,
                                         numericInput('Peff_caco2', label = NULL, value = 10)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Hepatic Clearance (L/h)')
                                  ),
                                  column(6,
                                         numericInput('Clh', label = NULL, value = 10)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Renal Clearance (L/h)')
                                  ),
                                  column(6,
                                         numericInput('Clr', label = NULL, value = 10)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Enterocyte Clearance (L/h)')
                                  ),
                                  column(6,
                                         numericInput('Clent', label = NULL, value = 0)
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Partition Coefficient Method')
                                  ),
                                  column(6,
                                         radioButtons("PCM", 
                                                            label = NULL, 
                                                            choices = list("Poulin & Theil" = "PT", 
                                                                           "Berezhkhovsky" = "bere"),
                                                            selected = "PT")
                                  )
                                ),
                                fluidRow(
                                  column(6,
                                         h5('Species')
                                  ),
                                  column(6,
                                         radioButtons("Species", 
                                                            label = NULL, 
                                                            choices = list("Human" = "human", 
                                                                           "Mouse" = "mouse",
                                                                           "Beagle" = "beagle",
                                                                           "Dog" = "dog"),
                                                            selected = "human")
                                  )
                                ),
                                p(' '),
                                actionButton("runButton", "Run simulation", width='250px'),
                                hr()
                              ),
                              mainPanel(
                                fluidRow(
                                ),
                                plotOutput("PK", height="1000px"),
                                br()
                              )
                            )
                   ),
                   tabPanel("Model description",
                            shinyUI(fluidPage(
                              withMathJax(),
                              h1('The Manchester PBPK model'),
                              p(),
                              br(),
                              p('The Manchester physiologically based pharmacokinetic (PBPK) model is designed as an educational tool. The aim is to provide a simple, open source, freely downloadable PBPK model for performing basic single subject bottom-up simulations in R.'),
                              em('This is a beta version of the software.'),
                              hr(),
                              h3('Basic equations of the PBPK model'),
                              p('This model is composed of two parts: a PBPK model describing the distribution, metabolism and elimination of the drug in the body and a compartmental absorption and transit (CAT) based model describing events following per oral drug administration in the gut.'),
                              br(),
                              
                              h4('PBPK model for drug distribution in the body'),
                              p('The PBPK model is a linear ordinary differential equations (ODE) model. This model is composed of 15 compartments, representing the lungs, brain, heart, kidneys, bones, muscles, stomach, spleen, liver, gut, pancreas, skin, fat, arterial and venous blood. PBPK model structure is shown in the image below, were red, blue and black-dashed arrows represent arterial and venous blood flow and clearances, respectively.'),
                              div(img(src = "PBPK_model.png", height = 400), style="text-align: center;"),
                              p('The drug distributes in all the organs. Clearance is supposed to happen only in liver and kidneys. The Equation describing drug distribution in all the tissues except the liver, lungs, arterial and venous blood is the following.'),
                              p(withMathJax('$$\\frac{dc_i}{dt} = \\frac{Q_i}{V_i} \\biggl(c_{art} - \\frac{c_i}{K_{i,b}/BP}  \\biggr)$$')),
                              p('\\(c_i\\) and \\(c_{art}\\) are the i-th organ and arterial blood concentration; \\(Q_i\\) and \\(V_i\\) are the i-th organ blood flow and volume; \\(K_{i,b}\\) is the tissue to blood partition coefficient and \\(BP\\) is the blood to plasma ratio.'),
                              p('Equations for lungs (l), liver (liv), kidneys (kid), arterial blood (art) and venous blood (ven) are the following.'),
                              p(withMathJax('$$\\frac{dc_l}{dt} = \\frac{Q_{all}}{V_l} \\biggl(c_{ven} - \\frac{c_l}{K_{l,b}/BP}  \\biggr)$$')),
                              p(withMathJax('$$\\frac{dc_{liv}}{dt} = \\frac{Q_{liv}}{V_{liv}} \\biggl(c_{art} - \\frac{c_{liv}}{K_{liv,b}/BP}  \\biggr) + \\frac{1}{V_{liv}}\\sum_{j\\in S} \\biggl( Q_j \\frac{c_j}{K_{j,b}/BP} \\biggr) - CL_h \\cdot \\frac{c_{liv}}{V_{liv}} + input_{GI}$$')),
                              p(withMathJax('$$\\frac{dc_{kid}}{dt} = \\frac{Q_{kid}}{V_{kid}} \\biggl(c_{art} - \\frac{c_{kid}}{K_{kid,b}/BP}  \\biggr) - CL_r \\cdot \\frac{c_{kid}}{V_{kid}}$$')),
                              p(withMathJax('$$\\frac{dc_{art}}{dt} = \\frac{Q_{all}}{V_{art}} \\biggl( \\frac{c_l}{K_{l,b}/BP} - c_{art}  \\biggr)$$')),
                              p(withMathJax('$$\\frac{dc_{ven}}{dt} = \\frac{1}{V_{ven}}\\sum_{j\\in D} \\biggl( Q_j \\frac{c_j}{K_{j,b}/BP} \\biggr) - \\frac{Q_{all}}{V_{ven}} c_{ven} $$')),
                              p('\\(S\\) is the set of the splanchnic organs (stomach, spleen, gut, pancreas); \\(D\\) is the set of the organs whose venous output enters directly the venous blood compartment (brain, heart, kidneys, bones, muscles, liver, skin, fat); \\(Q_{all}\\) is the cardiac output; \\(CL_h\\) and \\(CL_r\\) are the hepatic and renal clearances; \\(input_{GI}\\) is the input from the CAT based model.'),
                              br(),
                              
                              h4('Compartmental absorption and transit based model'),
                              p('Similarly to the PBPK model for drug distribution, the CAT based model is a compartmental model composed of a linear ODE system. The model represents the dissolution and transit out of the stomach and dissolution, transit and absorption occuring in the small intestine. In this model, the small intestine is divided in 6 different sections: one for the duodenum, two for the jejunum and three for the ileum. CAT based model structure is shown in the image below, where St, D, J, I and LI stand for stomach, duodenum, jejunum, ileum and large intestine and \\(M_{0}\\) is the drug dose.'),
                              div(img(src = "CAT_based_model.png", height = 300), style="text-align: center;"),
                              p('Equation for describing drug transit and dissolution in the stomach are reported below.'),
                              p(withMathJax('$$ \\frac{dx_{st,s}}{dt} = -k_{t,0} x_{st,s} - K_{st} x_{st,s} $$')),
                              p(withMathJax('$$ \\frac{dx_{st,d}}{dt} = -k_{t,0} x_{st,d} + K_{st} x_{st,s} $$')),
                              p('\\(x_{st,s}\\) and \\(x_{st,d}\\) are the solid and dissolved drug amount in the stomach; \\(k_{t,0}\\) is the time constant for drug output from the stomach and is calculated as the inverse of the gastric emptying time; \\(K_{st}\\) is the drug dissolution rate in the stomach and is described with the Noyes-Whitney model, as shown below.'),
                              p(withMathJax('$$ K_{st} = \\frac{3D}{\\rho h r} \\biggl( C_{st} - \\frac{x_{st,d}}{V_{st}} \\biggr) $$')),
                              p('\\(\\rho\\) is the density of the drug particle; \\(r\\) is the particle radius of the formulation; \\(h\\) is the effective thickness of the hydrodynamic diffusion layer and is calculated from \\(r\\) with the Hintz and Johnson model: \\(h=r\\) if \\(r<30 \\mu m\\), otherwise \\(h=30\\mu m\\). \\(C_{st}\\) is the drug solubility in the stomach and is calculated from the intrinsic drug solubility (\\(C_{int}\\)) as \\(C_{st}=C_{int} \\cdot \\alpha_{st} \\), where \\(\\alpha_{st}\\) is defined for neutral, monoprotic acidic and basic compounds using the Henderson Hasselbalch equation as follows.'),                              
                              p(withMathJax('$$ \\alpha_{neutral} = 1 $$')),
                              p(withMathJax('$$ \\alpha_{acid} = 1 + 10^{pH-pKa} $$')),
                              p(withMathJax('$$ \\alpha_{base} = 1 +10^{pKa-pH}$$')),
                              p('\\(pKa\\) is the drug dissociation constant while pH is the pH of the solvent (given stomach or intestine section where the drug dissolves).'),
                              p('In the Noyes-Whitney model, \\(D\\) is the drug diffusion coefficient and can be calculated from the Stokes-Einstein equation, as follows.'),
                              p(withMathJax('$$ D = \\frac{k_b T}{6\\pi\\eta_w R_h} $$')),
                              p('\\(k_b\\) is the Boltzmann constant, \\(T\\) is the absolute temperature of the body in Kelvin, \\(\\eta_w\\) is the viscosity of water at body temperature and \\(R_h\\) is the hydrodynamic radius of the diffusing drug. \\(R_h\\) is calculated as follows, assuming the drug molecule is spherical in shape.'),
                              p(withMathJax('$$ R_h = \\sqrt[3]{\\frac{3 mw}{4\\pi N_A \\rho}} $$')),
                              p('\\(mw\\) is the compound molecular weight while \\(N_A\\) is the Avogadro\'s number.'),
                              p('Equations for describing the transit, dissolution and absorption happening in the i-th section of the small interstine are reported below.'),
                              p(withMathJax('$$ \\frac{dx_{i,s}}{dt} = k_{t,i-1} x_{i-1,s} - k_{t,i} x_{i,s} - K_{i} x_{i,s} $$')),
                              p(withMathJax('$$ \\frac{dx_{i,d}}{dt} = k_{t,i-1} x_{i-1,d} - k_{t,i} x_{i,d} + K_{i} x_{i,s} - k_{a,i} x_{i,d}$$')),
                              p(withMathJax('$$ \\frac{dx_{i,ent}}{dt} = k_{a,i} x_{i,d} - CL_{ent} \\frac{x_{i,ent}}{V_{i,ent}} - Q_{i,ent} \\frac{x_{i,ent}}{V_{i,ent}}$$')),
                              p('\\(x_{i,s}\\), \\(x_{i,d}\\) and \\(x_{i,ent}\\) are the amount of solid and dissolved drug in the i-th compartment of the small interstine and the amount of drug in the i-th enterocytic compartment; \\(Q_{i,ent}\\) and \\(V_{i,ent}\\) are the blood flow and volume of the i-th section of the enterocytes; \\(K_{i}\\) is the drug dissolution rate in i-th segment of the small intestine and its definition is equivalent to \\(K_{st}\\); \\(CL_{ent}\\) is the clearance happening in the enterocytes compartments and is supposed to be equal for all the six sections. \\(k_{t,i}\\) is the transit time constant for the i-th small intestine compartment and is calculated as \\(k_{t,i} = (SITT\\cdot l_i/l_{tot})^{-1}\\), where \\(SITT\\) is the small intestinal transit time, \\(l_i\\) is the small intestine segment length and \\(l_{tot}\\) is the total length of small intestine. \\(k_{a,i}\\) is the absorption constant of the i-th compartment of the small intestine and is calculated from the effective jejunal permeability (\\(P_{eff}\\)) as \\(k_{a,i}=2 P_{eff} /R_{i}\\), where \\(R_{i}\\) is the radius of the intestinal compartment. \\(input_{GI}\\) in the PBPK equations is defined as follows.'),                            
                              p(withMathJax('$$ input_{GI} = \\sum_{i=1}^{6} Q_{i,ent} \\frac{x_{i,ent}}{V_{i,ent}}$$')),
                              hr(),
                              
                              h3('PBPK input parameters'),
                              p('Starting point is to provide the PBPK input parameters. Currently, the units of those parameters are considered to be standard and not customizable.'),
                              tags$div(
                                tags$ul(
                                  tags$li("mw [g/mol]: molecular weight"),
                                  tags$li("Type: select if neutral, monoprotic acid or base"),
                                  tags$li("logPow: octanol to water partition coefficient, used to calculate the partition coefficients"),
                                  tags$li("fup: fraction unbound in plasma, currently used to calculate the partition coefficient"),
                                  tags$li("BP: blood to plasma ratio, used to derive the plasma concentration from blood concentration"),
                                  tags$li("pKa: dissociation constant, used to calculate water solubility in gut sections"),
                                  tags$li("r [\\(\\mu\\)m]: radius of the particle size of the formulation, used to calculate the dissolution coefficient"),
                                  tags$li("Density of the formulation [g/L]: used to calculate the dissolution coefficient"),
                                  tags$li("Intrinsic solubility [mg/L]: used to calculate the dissolution coefficient"),
                                  tags$li("Peff [\\(10^{-4}\\)cm/s]: effective permeability across the gut layer"),
                                  tags$li("Hepatic clearance [L/h]: clearance referred to total concentration of drug into the liver"),
                                  tags$li("Renal clearance [L/h]: clearance referred to total concentration of drug into the kidneys"),
                                  tags$li("Enterocyte clearance [L/h]: clearance referred to total concentration of drug into the enterocytes, assumed to be equal for all the enterocytes sections"),
                                  tags$li("Partition coefficient method: select either the Poulin & Theil or Berezhkhovsky method"),
                                  tags$li("Species: select either human, mouse, beagle or dog"),
                                )
                              ),
                              p('In addition, the user can select the desired schedule.'),
                              tags$div(
                                tags$ul(
                                  tags$li("Route: currently only intra-venous (iv) and per oral (po) boluses are supported; for po route, drug can be administered both in solid and dissolved form."),
                                  tags$li("Schedule: currently once and twice daily dose administrations are supported"),
                                  tags$li("Dose [mg]"),
                                  tags$li("Days: for how many days the schedule needs to be repeated.")
                                )
                              ),
                              p('All the remaining physiological parameters are fixed to mean values. Values for these parameters and relative references can be found with the source code on github, in the data directory.'),
                              hr(),
                              
                              h3('References'),
                              p("The model code can be found at the following address: PUT GITHUB LINK"),
                              p("Useful references for the PBPK modelling:"),
                              tags$div(
                                tags$ul(
                                  tags$li(tags$a(href="https://doi.org/10.1038/psp.2013.41", "Jones, Rowland-Yeo 2013"), ": tutorial showing basic concepts of single subject and population PBPK models and in vitro to in vivo extrapolation."),
                                  tags$li(tags$a(href="https://doi.org/10.1002/jps.21798", "Berezhkovskiy 2009"), ": well explained theoretical bases of PBPK models."),
                                  tags$li(tags$a(href="https://doi.org/10.1016/S0378-5173(99)00147-7", "Yu Amidon 1999"), ": explanation of the compartmental absorption and transit model."),
                                  tags$li(tags$a(href="https://doi.org/10.1016/j.xphs.2018.10.033", "Grimstein et al. 2019"), ": FDA review on regulatory use of PBPK models."),
                                  tags$li(tags$a(href="https://doi.org/10.1007/s10928-018-9615-8", "Melillo et al. 2018"), ": publication from University of Manchester and University of Pavia groups considered as a reference for the equations of the CAT based model."),
                                )
                              ),
                              p("Regulatory guidances on the use of PBPK models can be found here."),
                              tags$div(
                                tags$ul(
                                  tags$li(tags$a(href="https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-reporting-physiologically-based-pharmacokinetic-pbpk-modelling-simulation_en.pdf", "European Medicines Agency (EMA), Committee for Medicinal Products for Human Use (CHMP), 2018")),
                                  tags$li(tags$a(href="https://www.fda.gov/media/101469/download", "Food and Drug Administration (FDA), Center for Drug Evaluation and Research (CDER), 2018")),
                                  tags$li(tags$a(href="https://www.oecd.org/chemicalsafety/risk-assessment/guidance-document-on-the-characterisation-validation-and-reporting-of-physiologically-based-kinetic-models-for-regulatory-purposes.pdf", "Organisation for Economic Co-operation and Development (OECD), 2021")),
                                  tags$li(tags$a(href="https://www.who.int/ipcs/methods/harmonization/areas/pbpk_models.pdf?ua=1", "International Programme on Chemical Safety (IPCS), World Health Organization (WHO), 2010")),
                                  )
                              ),
                              hr(),
                              
                              h3('Authors'),
                              tags$div(
                                tags$ul(
                                  tags$li("Hitesh Mistry, University of Manchester: conceptualization of the work and development of Shiny R app."),
                                  tags$li("Nicola Melillo, University of Manchester: development of PBPK software and Shiny R app."),
                                  )
                              ),
                              br(),
                              br(),
                              
                            )),
                            
                            
                            
                   )
)
)



## to insert latex equations
# p('Inline equation! \\(x^2\\)'),
## equation 
# p(withMathJax('$$x^2$$')), https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-reporting-physiologically-based-pharmacokinetic-pbpk-modelling-simulation_en.pdf


