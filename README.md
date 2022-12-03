# rnaseq_preprocess

![CI](https://github.com/ATpoint/sc_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.6-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000&logo=data%3Aimage%2Fjpeg%3Bbase64%2C%2F9j%2F4AAQSkZJRgABAQABLAEsAAD%2F2wBDABwcHBwcHDAcHDBEMDAwRFxEREREXHRcXFxcXHSMdHR0dHR0jIyMjIyMjIyoqKioqKjExMTExNzc3Nzc3Nzc3Nz%2F2wBDASIkJDg0OGA0NGDmnICc5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ub%2FwAARCAH0Ae8DASIAAhEBAxEB%2F8QAGgABAAIDAQAAAAAAAAAAAAAAAAQFAgMGAf%2FEADcQAQACAQIDBQYFAwQDAQAAAAABAgMEERIhMQUyQVFxEyJSYYGRM0KhsdEUcsEVI0NiNOHw8f%2FEABkBAQADAQEAAAAAAAAAAAAAAAABAgMEBf%2FEACIRAQEAAgIDAQEBAQEBAAAAAAABAhEDMRIhURNBYUIiUv%2FaAAwDAQACEQMRAD8A6QAAAAAAAAAAAAAAAAAAAAAAAAOnVqtnw163g0i3TaIltbgjpvP0ap7Qr%2BWk%2FWVvC%2FFbyY%2FVgKue0L%2BFYYTrs3lC351X9cVuKrHrM1slazttMxHRaq5Y2dr45zLoELJqb0vNYiNoeRq7eNYZfpE7ThDjV18ay2RqsU9d4T5z6bSBhGXHbpaGayQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJmIjeeSHk1uKnKnvT%2BiZLekXKTtMa75cePv2iFRk1ebJy34Y%2BSM1nF9Y5c3xa319I5UrM%2BvJFvrM9uk8PoiC8wkY3kyv9ZWve%2FemZ9WILqAAAANmH8Wn90Ogc%2Fh%2FFp%2FdDoGHL26eDqqrP8Ai29Wptz%2FAItvVqedl3WlAEIGVb3r3ZmGICRXVZI67SkV1VJ70bK8XnJYna4ret%2BdZ3ZKWJmJ3jkkU1OSve96Gk5Z%2FU7WQ0U1GO%2FLfafm3tZZelgBIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAj5tTjw8p528oTJvpFsnupEzERvKDl1tKcsfvT5%2BCBl1GTNPvTy8oaG2PH9c%2BXN%2F8tuTNkyzved%2Fl4NQNdMLd9gAAAAAAAAANmH8Wn90Ogc%2Fh%2FFp%2FdDoGHL26eDqqrP8Ai29Wptz%2FAItvVqedl3WlAEIAAAAAAG3HnyY%2Bk7x5S1CZbOhZ49RTJy6T5N6lSMeovTlbnDbHl%2BrSrIYUyVyRvWWbZYAAAAAAAAAAAAAAAAAAAAAAAAAAeTMVjeZ2iGGTLTFXivKnz6i%2BaefKvhC%2BOFrPPkmKRn1sz7uHlHmr%2BvOQdExk6cuWVy7AEqgAAAAAAAA2VxZb92sy3V0eefDb1lFyi0xt6iKJ8aDJ%2Ba0Qzjs%2Fzv8Aor54rfll8QcP4tP7odAg00NKWi3FPKd05lyZS9N%2BLG4z2qs%2F4tvVqWGTTcdptFtt2mdJk8JhxZYXa2kUbp0%2BWPDdrmlq96JhS42DEBCAAAAAAAAHtbWrO9Z2lYYtRF%2FdvylXC2OdiZV0K%2FDqJr7uTnHmsImJjeHTjlL0tKALJAAAAAAAAAAAAAAAAAAAAGjPnpgrvPOZ6Q81Gorgr52npClve2S02tO8y0ww37rLk5PH1GWTJfLbivLWDocluwAAAAAAZ0x3yTtSN1hi0ERzyzv8oVuUna%2BOFy6VsRNp2rG8pdNFmvzt7sfNbUpTHG1IiGTK8t%2FjbHhn9QqaHFXvzNv0Sa4sdO7WIbBncre2sxk6gAhYAAAAAAABrthx261hotpKz3J29UsVuMvaNKu%2BDLTw3j5NK6a74seTvR9Wd4viNKkSsmltHOnOP1RZiYnaeTG42doAEIAAAAG7FmtinbrXyaRMuvcFxW1bxxVneGSpxZbYp3jp4wtKXrkrxVdOGfkvKyAXSAAAAAAAAAAAAAAAANGfPXBXfrM9IZZctcNOO30hR5Mlst5vbrLTDDbLk5PH1Hl72vabWneZYg6HIAAAAAkYNPfNO8cq%2BaLddpkt9RorW1p4axvKxw6H82b7Qm4sOPDG1I%2Bvi2scuTfTpw4pPdeVrWscNY2h6DJsAAAAAAAAAAAAAAAAAAML46ZI2tDMBXZNNanOvOP1Rl0j5dPXJzjlLHLi%2BK2K0ZXpak8No2YsFQAAABsxZbYrbx08Yawl0Lil63rxV6MlVhyzit8p6wtImLRvHSXVhn5ReV6AukAAAAAAAAAAAAeWtFKza3KIeqjWZ%2FaW9nXux%2BsrY47qmeXjNtGfNbNfinp4Q0g6ZNOO3fugCUAAALTTaTh2yZY5%2BEK5ZSLY4XK6jVp9HNtr5eUeEea0iIiNo5Q9HPlla7McJj0AKrAAAAAAAAAAAAAAAAAAAAAAAAMb0rkjhtCty4bYp8481o8mItG084UywlRYphIzYJx%2B9XnX9kdzWWeqqAIQAAJOnzcE8Nu7P6IwmXV3EroRNNl4o9nbrHRLdeN3NrgCQAAAAAAAABje8Y6ze3SARdXn9nTgr3rfpCnZ5Mlsl5vbxYOrHHUcWeflQBZQAABZaTTdMuSPSP8AKuWWptbHG5XUZ6XS8G2TJHPwjyTwc1u%2FddmOMk1ABCwAAAAAAMbWrSN7zER80LJ2npMf5uKf%2BvM0J4or9tR%2Fx4%2FvKLftbVW7vDX0j%2BVvGo26cchbX6y3XJP05fs0Tnz263tP1lPibdtvs84q%2BcOGmZnrLw8Ebd1xV84OKvm4UPA27scNF717szHo2xqdRXpktH1k8E7doORr2hrK9Mk%2FXaUiva%2Bqr3orb6I8abdMKKnbUf8AJj%2B0pdO1tJbvTNfWP4Rqm1kI9NXpsncyVn6pHXohIAAABMRMbSrs%2BD2c8Ve7%2ByxJiJjaVcsdxFilG%2FPhnHO8d2Why2auqqAIQAA9iZrO8dYWuLJGSnF4%2BKpbsGT2d%2BfSerTjy1UyrQB0rgAAAAAAACr12be3sq9I6rDLkjFjm8%2BCgmZtM2nrLXjx%2FrDmy1NPAG7mAAAbcWK2W8Uj6hJv036TT%2B1tx27sfquGNKVpWKV6QycuWW67cMPGaAFVwAAAAY3vTHWb3naI8Zc%2Fq%2B1b33pp%2Fdr8XjP8Jk2bXGo1mDTfiW5%2BUdVLn7XzX5YY4I8%2BsqmZmZ3nnLxeYxXbO%2BTJlniyWm0%2FNgCyAAAAAAAAAAAAAABspmy4%2FwAO819JawFni7V1WPvTF4%2BcfwscXbGC3LLWaT94c2I8YnbtcWow5vwrxb0bnIaDF7bVUr4RO8%2FR17OzSYAIS8tWL1mtukqnJjnHbhlbtObF7SnLrHRnnjuIsVYdOUjmUAAAAWOmycdOGetUlU4r%2BzvFvutuvOHTx5bi8oA0SAAAAAxvaKVm89IjcFZrsu9oxR0jnPqgMrWm9ptPWWLrxmppw5ZbuwBKoAAutLh9lj3nvW6oOjw%2B0ycU9K%2FuuGPJl%2FHRw4f9UAYugAAAAas2bHgxzkyTtENlrRWJtadojnLktbq7arLvHcr3Y%2FymTaLWOr1mTV33tyrHSqIDVUAAAAAAAAAAAAAAAAAAAAB7ETM7QC%2F7Hw7Uvnnx5R9F206bFGDBTF5Rz9fFuZWrwAQAAIOqxbT7Sv1Q1zMRaJrPSVTkpOO81nwc%2FJjr2rYwAZKgACx01%2BKnDPWqubsF%2BDJHlPJfDLVTFoA6lwAAABB12Thxxjj837QnKTV3480%2BUcl%2BObrLly1ijAOlyAAB15QJejx8eXinpXmi3U2nGbulngxeyxxXx8fVuByW7d0mvQAJAAAeTMRG89IBTdr6nhrGmpPO3O3o55u1GWc%2Ba2WfzTy9PBpayaUoAkAAAAAAAAAAAAAAAAAAAAE%2Fs7D7bVV36V96fogOj7Hw8OG2aet52j0hGV9EXADJcAAAARdVj3rxx1j9kp5MRMbT4oym5oUwyvWaWms%2BDFxswAAAFtivx44s2IWkt1p9U1143c2vABZIADHJbgpN58I3c7MzM7yuNbfhw8PxTsp2%2FFPW3NzX3oAasAABdaPHwYYmetuaox09petPOXQxG0bQy5b%2FABvw4%2F0AYOkAAAAQu0MnstJeY6z7v3TVL2zfbHjx%2BczP2%2F8A1M7K54BqoAAAAAAAAAAMqUtktFKRvM9Ihi6LsjT1rinUT3rTtHpCLdEacPY8zG%2Be%2B3yr%2FKXHZGl87ff%2FANLQZ%2BVW0qbdj6ee7a0fb%2BETJ2NkjnivE%2Bsbfy6EPKmnG5tJqMH4lJiPPrCM7vrylWansvBm3tj9y3y6fZaZfUacuN%2BfTZtNbhyxt5T4S0LoAAe1rNrRWOsztDtsOOMOKuKPyxs5nsvD7XVRaelOf8OqUyq0AFEgAAAAAIOrpzi8ePJDWuavHjmFU5uWaqlAGaAAG3BbhyxP0WqlXFLcVYt5w34r%2FFsWQDZYABV6%2B29608o3%2B6vSdXbiz2%2BXJGdWE1I4uS7yoAsoAAnaGnFlm%2Fwx%2B62QtDXbFNvilNc2d3XZxTWIAo0AAAAHPdsz%2Fu46%2BVZl0Lm%2B2f8AyK%2F2%2FwCZWx7RVQA0VAAAAAAAAAAHS9kZq2wTh%2FNSf0lzTZiy3w3jJjnaYRZsjtxVaXtXFl2rm9y3n4StYmJjeGdi4AgAAYZMWPNSaZI3iXM63s%2B%2Bmnjp72Pz8vV1LyYi0TW0bxKZdIscKLTtDQTp59ri545%2FRW0pOS8Ur1tO0fVpKq6TsnDwaecs9bz%2BkLVhjpGPHXHXpWNmbOrgCAAAAAAAVGWvBkmvzW6v1ddrxbzhlyz1tFRQHOoAALLTW3xbeU7K1N0c96v1acd9piaA6VwGN54aWt5RIKDJPFe1vOZYA7Hn0AAB7EbzEeYL7BXhw0j5NpEbRsOOu%2BTQAJAAAAHP9tR7%2BO3nEw6BU9sY%2BLTxePy2%2FdOPaK5oBqqAAAAAAAAAAAAJWDWajTfh25eU84RQHSYO18N%2BWaOCfPrC1pemSvFSYtHnDhmzHlyYp4sdprPyVuKdu3HOYO18tOWasXjzjlK4wa7TZ%2BVLbT5TylSyp2lgIS8tWt6zS0bxPKYU2n7PnBruLrjrG9Z%2FwuhOwAQAAAAAAAACLq43pE%2BUpTTnjfFZXKblRVWA5FAABI0s7ZdvOEdtwztlr6rY9xMWoDrXGnUTtgv6NyNq5209vp%2B6ce1cuqpAHW4QABtwRxZqR84akjSRvqK%2F%2FeCMulse4vAHI7gAAAAABp1GL22C%2BL4o%2FVuAcLMTE7T1h4te1dN7LN7ase7f91U1lUAEgAAAAAAAAAAAAAAACbg1%2BpwcotxV8rc13p%2B1MGbauT%2Fbt8%2Bn3cuIuMTt3cTE84HH6fW59NO1J3r8M9HRaXX4dT7vdv5T%2FhS46TtOAVSAAAAAAAAMbxvS0fJkApQHEzAAGVJ2vE%2FNiAugjnA7Wgi6z%2Fx5%2BiU8mItG1o3j5pl1doym5pzg6H2WL4a%2FY9li%2BGv2bfr%2FAI5%2Fwv1zw6H2WL4a%2FY9li%2BGv2P1%2Fw%2FC%2FXPJei%2FHj0lbeyxfDX7PYx0rO9axE%2FKEXk3NLY8OrvbIBi3AAAAAAAAas%2BGmoxTiv0n9HH58N9PknFk6x%2BrtUXV6THq8fDblaOk%2BS2N0ixxw359Pl01%2BDLG3lPhLQ0VAAAAAAAAAAAAAAAAAAHvTnDwBd6LtS1ZjFqZ3jwt4x6r%2BJiY3jnEuFW%2FZ2unDaMGWfcnpPlP8ACmWPxMrpAFFgAAAAAAAFNbvT6vFx7Onwx9nns8fwx9mH5f6r4qgW%2Fs8fwx9j2eP4Y%2Bx%2BX%2BniqBb%2Bzx%2FDH2PZ4%2Fhj7H5f6eLKvdj0eg3WGGS00rxQza8sb45Bp%2FqJ8j%2BonyRgEn%2BonyP6ifJGASf6ifJlTNN7RWY6ojZjnbJHqCeAAAAAAAAAAADXlw489ODLG8Od1XZeXDvfD79f1h0wmXSNOEHXajQafU87Rw2%2BKFFqOzdRg51jjr5x1%2By8yRpXALIAAAAAAAAAAAAAAAAAAdP2XqpzYvZXn3qfrC0cdo8%2F9PqK5PDfafSXYs8otABVIAAAADG87UtPyBC%2Fq7fDB%2FV2%2BGEQcX65fXV%2BePxL%2Fq7fDB%2FV2%2BGEQP1y%2Bn54%2FEv%2Brt8MPY1dt%2B7CG9rztEfMnLl9PDH4uQHa5RjeN6zHyZAKwezG07PAAAHsTtMS8AWYxpO9In5MgAAAAAAAAGmdRjjPGnmdrTG8Nzj9XnnJq7Zaz0nl9EybRa7AQNBrI1WPa3K9esefzT0JAARc%2Bj0%2Bo%2FErz845Sps%2FZGWnvYJ448p5S6MTLYjTh70vjtw3iaz5Swdvkw4s1eHLWLR81PqOx4n3tPbb%2FrP8rzJGlANuXBlwW4ctZrLUsgAAAAAAAAAAAAAAdjosntdLjvPXbafpycc6Tse%2B%2BntX4bfurl0mLcBmsAAAANWedsVvRtRtVO2PbzlXO6xq2M9xXAOB1gADPFG%2BSsfOGDdp43yx8k4zdiMuloA9BxgAIGWNsktaRqI2tE%2BcI4AAAAJuCd6beTci6eec180oAAAAAAAAEbWZfYaa%2BTx22j1lxq%2F7Zy8qYI8fen%2FCgaYz0rWzFlvhyRkxztMOs0mrx6vHxV5WjrHk49sxZb4bxkxztMJs2Su3Fbo%2B0ceoiKZPdv8ApPosmWlgAAAGGTHTLSaZIi0S5TW6WdJl4Y51tzrLrlV2vSJ00WnrW0fqtjUVzIDRUAAAAAAAAAAAAX3Ys8ssen%2BVCvuxY%2FFn0%2Fyrl0mL0BmsAAAAIOrtzrX6pyqz24ssz9GPNdY6acU9tQDkdIAAl6SN7zbyhEWGkrtSbectOKbyU5L%2FAOUoB2uUABo1Eb1ifJDWN44qzCuAAAABnS3DeJWCsWGK3FSJBmAAAAAADRqsvsNPfJ4xHL18AcvrsvttVe3hE7R9EMGygAAtNL2pmw7Uy%2B%2FX9YVYWDssGr0%2Boj%2Fbtz8p5SkuE6c4TcXaGrxcovMx5TzUuKduuHPU7ZyxHv44n05fyzntqduWL9f%2FAEjxqdr5zvaurrkmNPjneKzvM%2FPyRs%2FaWpzxw78FZ8K%2Fyr0zFFoAugAAAAAAAAAAAAdJ2NXbBe3nb9oc267s%2FH7PSUjzjf7q5dJiaAzWAAAAYZLcFJt5KhN1d%2BlI9ZQnJzZbunTxTU2AMWgAAtsVeHHWvyVmOvHkivzW7o4J3WPLf4AOlgAAIGSvDeYT0bUV6W%2BgIoAAACRgttbh80d7EzE7x4Ash5W0WrFo8XoAAAACl7Zy7UphjxnefounI9oZfbaq8x0r7sfRbGe0VCAaKgAAAAAAAAAAAAAAAAAAAAAAANuDFObNXFH5p2drERWIiOkKPsjTddTaPlX%2FADK9Z5VaACqQAAEfU34Me0dbckZXU2mTd0gZL8d5swB59u%2FbskAAAAS9JXe038k9p09ODFHnPNud3HjrGOXO7oAuoAAMb14qzVkArXjdmrw338JaQAAAAScF9p4J%2BiUrYmYneFhS8XruDIAAAGnUZfY4L5fhjl6%2BDi5ned5dD2xl4cVcMfmnefSHOtMYrQBZAAAAAAAAAAAAAAAAAAAAAAAk6XTW1WWMdenjPlDDBgyajJGPHG8%2Fs63S6amlx%2Bzp18Z85Vt0mRupSuOkUpG0RG0MgZrAAAACqz5PaXmY6RyhM1OXhrwR1n9lc5ubP%2FmN%2BLH%2BgDnbAADPFTjvFWCdpKbROSfHlC%2FHju6VzuptMAdzkAAAAAAa8tOOnzhAWaFmpw25dJBpAAAAbcV%2BC3PpLUAsxGw5PyT9EkAAHKdpZZyau0eFfdhXu2y4MOaNstIt6qzL2PhtzxWmk%2BU84XmUVsc4LDL2Zq8XOK8cf9UG1bUna0TE%2FNfaGIAAAAAAAAAAAAAAAAAAJeDRanUdyu0ec8oBETtLoM2qni7tPin%2FAAuNN2VhxbWy%2B%2Fb9Fr05QpcviZGjBp8WmpwYo2858ZbwUWAAAAGN7xjrNpZdOasz5faW2jux0Z8mfjF8Md1qvab2m09ZYg4nUAAAAypWb2iseK3rWKxFY6Qi6XHtHtJ8eiW6%2BHDU25%2BTLd0ANmQAAAAAAwyU467ePgzAVvR4k56bTxx9UYAAAAHqbiyccbT1hBexMxO8Ashrx5IvHzbAAAGN8dMkbXrFo%2BcMgFfk7M0mTnFZrP8A1lCv2LH%2FAB5PvC9E7qNOYv2Tqq93ht6T%2FKNbQ6unXHP05%2Fs7AT5U04i2LLTv1mPWGt3bCceO3erE%2BsJ80acOOznSaWeuKv2hrnQaOf8AjhPmacgOsnszRT%2BT9ZY%2F6Xo%2Fhn7yeUNOVHVf6Vo%2Fhn7yf6Xo%2Fhn7yeUNOVHWx2boo%2F4%2F1lsjRaSvTFX7bnkacc20wZsncpafSHZ1xY6d2sR6QzR5mnK07L1l%2BtYr6ym4%2Bxo65cn0rH8r0R5VOkTDodLh50pEz5zzSwVSAAAAAAAiajPw%2B5Tr4yrllMZupxxtuow1Obf%2FAG6fVDBxZZXK7rrxx1NACqQABsxY5yXivh4tazwYvZ059Z6tOPDyqmeWo3xERG0AO1ygAAAAAAAAAExExtKBkpNLbJ7DJSL12BXj2YmJ2l4AAAADKJms7wm48kXj5oD2JmJ3gFkNOPLF%2BVuUtwAAAAAAAAAAAAAAAAAAAAAAAAAAAIWbU%2Flx%2FdXLKYzdWxxt6ZZ8%2FD7lOvjPkgA4s87ld1044yQAVWAAAbcOKcttvCOqZN3ULde27TYuKfaW6R0T3kRERtHSHruwx8ZpyZZbuwBZUAAAAAAAAAAABpy4%2BON46whrJHzYt%2Fer1BEAAAAAAScebblf7owCziYmN4ECmS1OiXTLW%2FLpINgAAAAAAAAAAAAAAAAAAAAADG960je07I%2BTU1rypznz8EG1rXne07yxz5ZPUa48dvbblz2yco5VaActytu63kk9QAQkAAB7Ws2nhrzmQe0pOS3DVa46Rjrw1Y4cUYq7eM9W12cfH4zdc2ee%2FUAGrMAAAAAAAAAAAAAAABHy4t%2Fer9kRZtGXDxe9XqCGPejwAAAAAAG6ma1eU84Sa5aX6TzQAFmIFct69Jb66is96NgSBjF626SyAAAAAAAAAAABqtmx06z9kWydpk302nTnKDfVz%2BSPujXyXv3p3ZZc0nTScVvadfU0ryr70oWTNfJ3p5eTWOfLkuTbHCQAUWAAAAAe1ra88NY3kCIm07RzmVlhwxijeedpe4cMYo8582518fHr3XPnnv1ABsyAAAAAAAAAAAAAAAAAAAAasmKL845Sh2rNZ2lYsb0reNpBXDZfHanzjzawAAAAAAAAGcZL16SwAb4z3jrtLZGojxhEATYz0%2BbL22PzQAFhGSk8olmr6d%2BPVYA1Wz4qTwzPOGudVjjpvKHn%2FFt6tTly5spbHROOa2mzrPhq1W1WWem0I4zvLlf6vMMfjK2S9u9MyxBS3fawAAAAAAAACVi00296%2FKE443K6iLlJ2048Vss7V6eayx4q442r92cVisbVjaHrsw45i5s87QBooAAAAAAAAAAAAAAAAAAAAAAAAdeqNfB40%2BySArZiYnaXixtSt496EW%2BG1edecA0D14AAAAAAAAAADOnfj1WCvp349VgCqz%2Fi29Wptz%2Fi29WpwZ912Y9QAVSAAAAAAA9rW1p2rG8g8Z0x3yTtWErHpfHJ9oTIiKxtWNob4cNvbLLlk6aMWnrj5zzlIB0zGT1GFtvYAlAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADC2Ot%2BsI1sFo7vNMAVsxMcpeLKa1tymN2m2nrPdnYEMbrYbx4b%2BjVMTHUHgAAAAAM6d%2BPVYK%2Bnfj1WAKrP8Ai29Wptz%2FAItvVqcGfddmPUAFUgAAyrS9%2B7EykV0t5707LTC3qIuUnaKzrjvfuxusKabHXrG8%2FNv225Q2x4PrK8vxCppPHJP0hLrWtI2rGzIb44THplcrewBZUAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAeTET15vQGucOOfDZrnTx4SkAIk6e%2FhMMJw5I8E4BX%2BzvHhLzht5SsQEClZ445eKeAKzNS85bbRM82MYMs%2FllajG8Mt3tr%2BtVkabLPhs2RpLfmtEJ4mcOKLy5IsaSkdZmW6uHHXpWGwXmEnUVuVv8AQBZUAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB%2F9k%3D)](https://sylabs.io/docs/)  

## Introduction

**rnaseq_preprocess** is a Nextflow pipeline for RNA-seq quantification with `salmon`. The processing steps are `fastqc` first, then quantification with `salmon`, aggregation to gene level with `tximport` and a small summary report with `MultiQC`. Multiple fastq files per sample are supported. These technical replicates will be merged prior to quantification. Optional trimming to a fixed read length is possible. The pipeline is containerized via Docker and Singularity. Outputs can be found in `rnaseq_preprocess_results/` including command lines and software versions. The expected Nextflow version is 21.10.6.

Run the test profile to see which output is being produced. Downloading the Docker image may take a minute or two:

```bash
NXF_VER=21.10.6 nextflow run atpoint/rnaseq_preprocess -r main -profile docker,test_with_existing_idx,test_resources
```

## Details

**Indexing**

The indexing step must be run first and separately using the `--only_idx` flag. For this we need a reference transcriptome (gzipped), a reference genome as decoy (gzipped) and a GTF annotation file (gzipped).

`--only_idx`: trigger the indexing process  
`--idx_name`: name of the produced index, default `idx`  
`--idx_dir`: name of the directory inside `rnaseq_preprocess_results/` storing the index, default `salmon_idx`  
`--idx_additional`: additional arguments to `salmon index` beyond the defaults which are `--no-version-check -t -d -i -p --gencode`  
`--txtome`: path to the gzipped transcriptome fasta  
`--genome`: path to the gzipped genome fasta  
`--gtf`: path to the gzipped GTF file  
`--transcript_id`: name of GTF column storing transcript ID, default `transcript_id`  
`--transcript_name`: name of GTF column storing transcript name, default `transcript_name`  
`--gene_id`: name of GTF column storing gene ID, default `gene_id`  
`--gene_name`: name of GTF column storing gene name, default `gene_name`  
`--gene_type`: name of GTF column storing gene biotype, default `gene_type`  

For the indexing process, 30GB of RAM and 6 CPUs are required/hardcoded. On our HPC we use:  

```bash
NXF_VER=21.10.6 nextflow run atpoint/rnaseq_preprocess -r main  -profile singularity,slurm --only_idx \
    --genome path/to/genome.fa.gz --txtome path/to/txtome.fa.gz --gtf path/to/foo.gtf.gz \
    -with-report indexing_report.html -with-trace indexing_report.trace -bg > indexing_report.log
```    

**Quantification/tximport**

The pipeline runs via a [samplesheet](./test/samplesheet.csv) which is a CSV file with the columns:
`sample,r1,r2,libtype`. The first column is the name of the sample, followed by the paths to the R1 and
R2 files and the salmon [libtype](https://salmon.readthedocs.io/en/latest/library_type.html). If R2 is left blank
then single-end mode is triggered for that sample. Multiple fastq files (lane/technical replicates) are supported.
These must have the same sample column and will then be merged prior to quantification. Optionally, a `seqtk` module can
trim reads to a fixed read length, triggered by `--trim_reads` with a default of 75bp, controlled by `--trim_length`. 
The quantification then runs with the salmon options `--gcBias --seqBias --posBias` (for single-end without `--gcBias`). 
Transcript abundance estimates from `salmon` are then summarized to the gene level using [tximport](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#Salmon) with its `lengthScaledTPM` option. That means returned gene-level counts are already corrected for average transcript length and can go into any downstream DEG analysis, for example with `limma`. Both a matrix of counts and effective gene lengths is returned.

Other options:

`--idx`: path to the salmon index folder  
`--tx2gene`: path to the tx2gene map matching transcripts to genes  
`--samplesheet`: path to the input samplesheet  
`--trim_reads`: logical, whether to trim reads to a fixed length  
`--trim_length`: numeric, length for trimming  
`--quant_additional`: additional options to `salmon quant` beyond `--gcBias --seqBias --posBias`  

We hardcoded 30GB RAM and 6 CPUs for the quantification. On our HPC we use:

```bash
NXF_VER=21.10.6 nextflow run atpoint/rnaseq_preprocess -r main -profile singularity,slurm \
    --idx path/to/idx --tx2gene path/to/tx2gene.txt --samplesheet path/to/samplesheet.csv \
    -with-report quant_report.html -with-trace quant_report.trace -bg > quant_report.log
```

**Other options**

`--merge_keep`: logical, whether to keep the merged fastq files  
`--merge_dir`: folder inside the output directory to store the merged fastq files  
`--trim_keep`: logical, whether to keep the trimmed fastq files  
`--trim_dir`: folder inside the output directory to store the trimmed fastq files  
`--skip_fastqc`: logical, whether to skip `fastqc`  
`--only_fastqc`: logical, whether to only run `fastqc` and skip quantification  
`--skip_multiqc`: logical, whether to skip `multiqc`  
`--skip_tximport`: logical, whether to skip the `tximport` process downstream of the quantification  
`--fastqc_dir`: folder inside the output directory to store the fastqc results  
`--multiqc_dir`: folder inside the output directory to store the multiqc results  
