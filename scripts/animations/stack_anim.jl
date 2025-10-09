
# function stack_gifs_vertical(gif1::String, gif2::String, out::String; width=400)
#     # Commande ffmpeg : on redimensionne les deux gifs à la même largeur
#     cmd = `ffmpeg -y -i $gif1 -i $gif2 -filter_complex \
#         "[0:v]scale=$width:-1[v0]; \
#          [1:v]scale=$width:-1[v1]; \
#          [v0][v1]vstack=2[v]" \
#         -map "[v]" $out`
#     run(cmd)
# end

# function stack_gifs_vertical(gif1::String, gif2::String, out::String; width=400)
#     # Fichiers temporaires
#     palette = tempname() * ".png"

#     # 1. Générer une palette optimisée
#     run(`ffmpeg -i $gif1 -i $gif2 -filter_complex \
#         "[0:v]scale=$width:-1[v0]; \
#          [1:v]scale=$width:-1[v1]; \
#          [v0][v1]vstack=2[v]; \
#          [v]palettegen=stats_mode=diff" \
#         -y $palette`)

#     # 2. Créer le gif final en utilisant la palette
#     run(`ffmpeg -i $gif1 -i $gif2 -i $palette -filter_complex \
#         "[0:v]scale=$width:-1[v0]; \
#          [1:v]scale=$width:-1[v1]; \
#          [v0][v1]vstack=2[v]; \
#          [v][2:v]paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle" \
#         -y $out`)
# end

# function stack_gifs_vertical(gif1::String, gif2::String, out::String; width=400)
#     palette = tempname() * ".png"

#     # 1. Générer une palette avec filtres nets
#     run(`ffmpeg -i $gif1 -i $gif2 -filter_complex \
#         "[0:v]scale=$width:-1:flags=lanczos[v0]; \
#          [1:v]scale=$width:-1:flags=lanczos[v1]; \
#          [v0][v1]vstack=2[v]; \
#          [v]palettegen=stats_mode=diff" \
#         -y $palette`)

#     # 2. Créer le gif final en appliquant la palette
#     run(`ffmpeg -i $gif1 -i $gif2 -i $palette -filter_complex \
#         "[0:v]scale=$width:-1:flags=lanczos[v0]; \
#          [1:v]scale=$width:-1:flags=lanczos[v1]; \
#          [v0][v1]vstack=2[v]; \
#          [v][2:v]paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle" \
#         -y $out`)
# end

function stack_gifs_vertical_noscale(gif1::String, gif2::String, out::String)
    palette = tempname() * ".png"

    # 1. Générer la palette
    run(`ffmpeg -i $gif1 -i $gif2 -filter_complex \
        "[0:v][1:v]vstack=2[v]; \
         [v]palettegen=stats_mode=diff" \
        -y $palette`)

    # 2. Appliquer la palette pour obtenir un GIF net
    run(`ffmpeg -i $gif1 -i $gif2 -i $palette -filter_complex \
        "[0:v][1:v]vstack=2[v]; \
         [v][2:v]paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle" \
        -y $out`)
end




# Exemple d’utilisation :
#stack_gifs_vertical_noscale("images/GTT/anim_test1.gif", "images/GTT/anim_test3.gif", "images/GTT/stacked_test.gif")

#stack_gifs_vertical_noscale("images/GTT/anim_MUSCL_u.gif", "images/GTT/anim_MUSCL_D.gif", "images/GTT/stacked_MUSCL.gif")

stack_gifs_vertical_noscale("mMUSCL_u.mp4", "mMUSCL_D.mp4", "anim_MUSCL.mp4")