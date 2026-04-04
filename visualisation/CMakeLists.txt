SET(DGTAL_TOOLS_CONTRIB
   displayTgtCoverAlphaTS
   displayLineSegments
   displaySet2dPts
)

SET(DGTAL_TOOLS_CONTRIB_POLY
   graphViewer
   polyMeshEdit
   polyMeshColorize
)





FOREACH(FILE ${DGTAL_TOOLS_CONTRIB})
  add_executable(${FILE} ${FILE}.cpp)
  target_link_libraries (${FILE} DGtal::DGtal ${DGtalToolsContribLibDependencies})
  install(TARGETS ${FILE} 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
ENDFOREACH(FILE)


FOREACH(FILE ${DGTAL_TOOLS_CONTRIB_POLY})
  add_executable(${FILE} ${FILE}.cpp)
  target_link_libraries (${FILE}  polyscope DGtal::DGtal  ${DGtalToolsContribLibDependencies})
  install(TARGETS ${FILE} 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
ENDFOREACH(FILE)


